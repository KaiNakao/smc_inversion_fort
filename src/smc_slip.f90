module smc_slip
    use omp_lib
    implicit none
contains
    subroutine slip_BoxMuller(p, q)
        implicit none
        double precision, intent(inout) :: p, q
        double precision :: pi, r, s
        pi = 2d0*asin(1d0)
        !   uniform random numbers between 0 and 1
        call random_number(r)
        call random_number(s)
        !   Gaussian random numbers,
        !   with weights proportional to eˆ{-pˆ2/2}and eˆ{-qˆ2/2}
        p = sqrt(-2d0*log(r))*sin(2d0*pi*s)
        q = sqrt(-2d0*log(r))*cos(2d0*pi*s)
    end subroutine slip_BoxMuller

    double precision function cdf_norm(x, mu, sigma2)
        implicit none
        double precision, intent(in) :: x, mu, sigma2
        ! cumulative distribution function of normal distribution
        cdf_norm = (1d0 + erf((x - mu)/sqrt(2d0*sigma2)))/2d0
    end function cdf_norm

    double precision function pdf_norm(x, mu, sigma2)
        implicit none
        double precision, intent(in) :: x, mu, sigma2
        double precision :: pi
        pi = 2d0*asin(1d0)
        !  probability distribution function of normal distribution
        pdf_norm = exp(-((x - mu)**2d0)/(2d0*sigma2))/sqrt(2d0*pi*sigma2)
    end function pdf_norm

    double precision function slip_calc_likelihood(svec, dvec, &
                                                   sigma2_full, gmat, &
                                                   log_sigma_sar2, log_sigma_gnss2, &
                                                   nsar, ngnss, gsvec, &
                                                   nobs, ndof)
        implicit none
        double precision, intent(in) :: svec(:), dvec(:), sigma2_full(:), &
            gmat(:, :), log_sigma_sar2, log_sigma_gnss2
        integer, intent(in) :: nsar, ngnss, nobs, ndof
        double precision, intent(inout) :: gsvec(:)
        integer :: i

        do i = 1, nobs
            gsvec(i) = 0d0
        end do
        call dgemv('n', nobs, 2*ndof, 1d0, gmat, &
                   nobs, svec, 1, 0d0, gsvec, 1)

        slip_calc_likelihood = nsar*log_sigma_sar2/2d0 + 3d0*ngnss*log_sigma_gnss2/2d0
        do i = 1, nobs
            slip_calc_likelihood = slip_calc_likelihood + &
                                   (dvec(i) - gsvec(i))**2d0/(2d0*sigma2_full(i))
        end do
    end function slip_calc_likelihood

    double precision function slip_calc_prior(svec, alpha2_full, theta, nplane, nxi, neta, &
                                              lmat_index, lmat_val, lsvec, nnode)
        implicit none
        double precision, intent(in) :: svec(:), alpha2_full(:), &
            lmat_val(:, :), theta(:)
        integer, intent(in) :: lmat_index(:, :), nnode, nplane, nxi, neta
        double precision, intent(inout) :: lsvec(:)
        integer :: n, i, j
        double precision :: tmp, log_alpha2
        n = 2*nnode

        do i = 1, n
            lsvec(i) = 0d0
        end do
        do i = 1, n
            tmp = 0d0
            do j = 1, 5
                tmp = tmp + lmat_val(j, i)*svec(lmat_index(j, i))
            end do
            lsvec(i) = tmp
        end do

        slip_calc_prior = 0d0
        do i = 1, nplane
            log_alpha2 = theta(8*i)
            slip_calc_prior = slip_calc_prior + &
                              2d0*(nxi + 1)*(neta + 1)*log_alpha2/2d0
        end do
        do i = 1, n
            slip_calc_prior = slip_calc_prior + lsvec(i)**2d0/(2d0*alpha2_full(i))
        end do
    end function slip_calc_prior

    subroutine slip_gen_init_particles(particles, likelihood_ls, prior_ls, nparticle, &
                                       ndim, llmat, max_slip, dvec, &
                                       sigma2_full, alpha2_full, theta, nplane, nxi, neta, &
                                       gmat, log_sigma_sar2, &
                                       log_sigma_gnss2, nsar, ngnss, nobs, ndof, &
                                       lmat_index, lmat_val, nnode, particle_cur, gsvec, lsvec)
        implicit none
        double precision, intent(in) :: max_slip, dvec(:), &
            lmat_val(:, :), llmat(:, :), gmat(:, :), sigma2_full(:), &
            alpha2_full(:), log_sigma_sar2, log_sigma_gnss2, theta(:)
        integer, intent(in) :: nnode, ndof, nsar, ngnss, nobs, &
                               nparticle, ndim, lmat_index(:, :), nplane, nxi, neta
        double precision, intent(inout) ::  particles(:, :), &
            likelihood_ls(:), prior_ls(:), particle_cur(:), gsvec(:), lsvec(:)
        integer :: iparticle, idim, jdim, cnt
        double precision :: x, fx, f1, f0
        double precision :: y_i, y_prev, err
        double precision :: v1, v2  ! box muller
        double precision :: sigma2_i, mu_i

        do idim = 1, ndim
            particle_cur(idim) = 0d0
        end do
        !   Gibbs sampling from Truncated Multi Variate Normal distribution
        !   obtain nparticle samples
        do iparticle = 1, nparticle
            !  component - wise update
            do idim = 1, ndim
                !  mean and variance for conditional probability distribution
                sigma2_i = &
                    (llmat(idim, idim)/alpha2_full(idim))**(-1d0)
                mu_i = 0d0
                ! loop for (jdim != idim)
                do jdim = 1, idim - 1
                    mu_i = mu_i - llmat(idim, jdim)/alpha2_full(idim)* &
                           particle_cur(jdim)*sigma2_i
                end do
                do jdim = idim + 1, ndim
                    mu_i = mu_i - llmat(idim, jdim)/alpha2_full(idim)* &
                           particle_cur(jdim)*sigma2_i
                end do

                ! sampling x from conditional distribution
                call slip_BoxMuller(v1, v2)
                x = mu_i + v1*sqrt(sigma2_i)
                fx = cdf_norm(x, mu_i, sigma2_i)
                f1 = cdf_norm(max_slip, mu_i, sigma2_i)
                f0 = cdf_norm(0d0, mu_i, sigma2_i)

                ! solve F(y) = (F(1) - F(0)) F(x) + F(0) by Newton's method
                ! where F(:) is CDF of normal distribution
                y_i = mu_i
                err = 1d0
                y_prev = y_i + err
                cnt = 0
                do while (err > 1d-8 .and. cnt < 100)
                    y_i = y_i - (cdf_norm(y_i, mu_i, sigma2_i) - ((f1 - f0)*fx + f0))/ &
                          pdf_norm(y_i, mu_i, sigma2_i)
                    err = abs(y_i - y_prev)
                    y_prev = y_i
                    cnt = cnt + 1
                end do
                particle_cur(idim) = y_i
            end do
            do idim = 1, ndim
                particles(idim, iparticle) = particle_cur(idim)
            end do
        end do

!$omp parallel do private(iparticle, idim, particle_cur, gsvec, lsvec)
        do iparticle = 1, nparticle
            do idim = 1, ndim
                particle_cur(idim) = particles(idim, iparticle)
            end do
            !   calculate negative log likelihood and prior
            likelihood_ls(iparticle) = &
                slip_calc_likelihood(particle_cur, dvec, sigma2_full, &
                                     gmat, log_sigma_sar2, log_sigma_gnss2, &
                                     nsar, ngnss, gsvec, nobs, ndof)
            prior_ls(iparticle) = slip_calc_prior(particle_cur, alpha2_full, theta, nplane, nxi, neta, &
                                                  lmat_index, lmat_val, lsvec, nnode)
        end do
!$omp end parallel do
    end subroutine

    subroutine slip_calc_mean_std_vector(vec, size, mean, std)
        implicit none

        double precision, intent(in) :: vec(:)
        integer, intent(in) :: size
        double precision, intent(inout) :: mean, std
        integer :: i
        mean = 0d0
!$omp parallel do private(i) reduction(+:mean)
        do i = 1, size
            mean = mean + vec(i)
        end do
!$omp end parallel do
        mean = mean/size
        std = 0d0
!$omp parallel do private(i) reduction(+:std)
        do i = 1, size
            std = std + (vec(i) - mean)**2d0
        end do
!$omp end parallel do
        std = std/size
        std = sqrt(std)
    end subroutine slip_calc_mean_std_vector

    double precision function slip_find_next_gamma(gamma_prev, likelihood_ls, &
                                                   weights, neglog_evidence, &
                                                   nparticle)
        implicit none
        double precision, intent(in) :: gamma_prev, likelihood_ls(:)
        integer, intent(in) :: nparticle
        double precision, intent(inout) :: neglog_evidence, weights(:)
        integer :: iparticle
        double precision :: min_likelihood, cv_threshold, lower, upper, &
            err, gamma, diff_gamma, mean, std, cv, likelihood, evidence
        ! find minimum of negative log likelihood
        min_likelihood = likelihood_ls(1)
!$omp parallel do private(iparticle) reduction(min : min_likelihood)
        do iparticle = 1, nparticle
            min_likelihood = min(min_likelihood, likelihood_ls(iparticle))
        end do
!$omp end parallel do

        cv_threshold = 1d0
        ! binary search for the next gamma
        ! such that c.o.v of the weight is equivalent to cv_threashold
        lower = gamma_prev
        upper = 1d0
        err = 1d0
        gamma = 1d0
        do while (err > 10d-8)
            gamma = (lower + upper)/2d0
            diff_gamma = gamma - gamma_prev
!$omp parallel do private(iparticle, likelihood)
            do iparticle = 1, nparticle
                likelihood = likelihood_ls(iparticle)
                ! extract min_likelihood to avoid underflow of the weight
                weights(iparticle) = &
                    exp(-diff_gamma*(likelihood - min_likelihood))
            end do
!$omp end parallel do

            ! calculate c.o.v = mean/std
            call slip_calc_mean_std_vector(weights, nparticle, mean, std)
            cv = std/mean
            if (cv > cv_threshold) then
                upper = gamma
            else
                lower = gamma
            end if
            err = abs(cv - cv_threshold)
            if (abs(gamma - 1) < 1d-8) then
                exit
            end if
        end do

        ! Calulate S_j(mean of the weight)
        evidence = 0d0
        diff_gamma = gamma - gamma_prev
!$omp parallel do private(iparticle, likelihood) reduction(+ : evidence)
        do iparticle = 1, nparticle
            likelihood = likelihood_ls(iparticle)
            evidence = evidence + exp(-diff_gamma*(likelihood - min_likelihood))
        end do
!$omp end parallel do
        evidence = evidence/nparticle
        neglog_evidence = -log(evidence) + diff_gamma*min_likelihood

        slip_find_next_gamma = gamma
    end function slip_find_next_gamma

    subroutine slip_normalize_weights(weights, nparticle)
        implicit none
        double precision, intent(inout) :: weights(:)
        integer, intent(in) :: nparticle
        double precision :: sum
        integer :: iparticle
        sum = 0
!$omp parallel do private(iparticle) reduction(+ : sum)
        do iparticle = 1, nparticle
            sum = sum + weights(iparticle)
        end do
!$omp end parallel do
!$omp parallel do private(iparticle)
        do iparticle = 1, nparticle
            weights(iparticle) = weights(iparticle)/sum
        end do
!$omp end parallel do
    end subroutine

    subroutine slip_calc_mean_particles(particles, weights, &
                                        mean, nparticle, ndim)
        implicit none
        double precision, intent(in) :: particles(:, :), weights(:)
        integer, intent(in) :: nparticle, ndim
        double precision, intent(inout) :: mean(:)
        integer :: iparticle, idim
        double precision ::weight
        do idim = 1, ndim
            mean(idim) = 0d0
        end do
!$omp parallel do private(iparticle, weight, idim) reduction(+:mean)
        do iparticle = 1, nparticle
            weight = weights(iparticle)
            do idim = 1, ndim
                mean(idim) = mean(idim) + weight*particles(idim, iparticle)
            end do
        end do
!$omp end parallel do
    end subroutine slip_calc_mean_particles

    subroutine slip_calc_cov_particles(particles, weights, mean, &
                                       cov, nparticle, ndim)
        implicit none
        double precision, intent(in) :: particles(:, :), weights(:), mean(:)
        integer, intent(in) :: nparticle, ndim
        double precision, intent(inout) :: cov(:, :)
        integer :: iparticle, idim, jdim, ierr
        double precision weight, di, dj
        do jdim = 1, ndim
            do idim = jdim, ndim
                cov(idim, jdim) = 0d0
            end do
        end do

!$omp parallel do private(iparticle, weight, idim, di, jdim, dj) &
!$omp reduction(+:cov)
        do iparticle = 1, nparticle
            weight = weights(iparticle)
            do jdim = 1, ndim
                dj = (particles(jdim, iparticle) - mean(jdim))
                do idim = jdim, ndim
                    di = (particles(idim, iparticle) - mean(idim))
                    cov(idim, jdim) = cov(idim, jdim) + weight*di*dj
                end do
            end do
        end do
!$omp end parallel do
        do jdim = 1, ndim
            do idim = jdim, ndim
                cov(idim, jdim) = cov(idim, jdim)*0.04d0
            end do
        end do

        ! LAPACK function for LU decomposition of matrix
        call dpotrf('L', ndim, cov, ndim, ierr)
        return
    end subroutine slip_calc_cov_particles

    subroutine slip_resample_particles(nparticle, weights, assigned_num)
        implicit none

        double precision, intent(in) :: weights(:)
        integer, intent(in) :: nparticle
        integer, intent(inout) :: assigned_num(:)
        integer :: iparticle
        double precision :: d_nparticle, u
        d_nparticle = nparticle
        ! systematic residual resampling
        call random_number(u)
        u = u/d_nparticle
        do iparticle = 1, nparticle
            assigned_num(iparticle) = &
                floor((weights(iparticle) - u)*d_nparticle) + 1
            u = u + assigned_num(iparticle)/d_nparticle - weights(iparticle)
        end do
    end subroutine slip_resample_particles

    subroutine slip_mcmc_sampling(gamma, particles, particles_new, &
                                  likelihood_ls, likelihood_ls_new, prior_ls, &
                                  prior_ls_new, cov, assigned_num, id_start, &
                                  nparticle, ndim, dvec, sigma2_full, alpha2_full, &
                                  theta, nplane, nxi, neta, &
                                  gmat, log_sigma_sar2, log_sigma_gnss2, nsar, &
                                  ngnss, nobs, ndof, lmat_index, lmat_val, &
                                  nnode, max_slip, st_rand_ls, metropolis_ls, &
                                  particle_cur, particle_cand, st_rand, gsvec, lsvec)
        implicit none
        double precision, intent(in) :: max_slip, dvec(:), &
            lmat_val(:, :), gmat(:, :), sigma2_full(:), alpha2_full(:), theta(:), &
            log_sigma_sar2, log_sigma_gnss2, gamma, cov(:, :)
        integer, intent(in) :: nnode, ndof, nsar, ngnss, nobs, nplane, nxi, neta, &
                               nparticle, ndim, lmat_index(:, :), assigned_num(:)
        double precision, intent(inout) ::  particles(:, :), &
            particles_new(:, :), likelihood_ls(:), &
            prior_ls(:), &
            likelihood_ls_new(:), prior_ls_new(:), &
            st_rand_ls(:, :), metropolis_ls(:), particle_cur(:), &
            particle_cand(:), st_rand(:), gsvec(:), lsvec(:)

        integer, intent(inout) :: id_start(:)
        integer :: iparticle, jparticle, idim, nassigned, istart
        double precision ::likelihood_cur, likelihood_cand, &
            prior_cur, prior_cand, post_cur, post_cand, metropolis
        double precision :: v1, v2 ! box muller
        do iparticle = 1, nparticle
            do idim = 1, ndim
                call slip_BoxMuller(v1, v2)
                st_rand_ls(idim, iparticle) = v1
            end do
            call random_number(metropolis)
            metropolis_ls(iparticle) = metropolis
        end do

        id_start(1) = 1
        do iparticle = 1, nparticle - 1
            id_start(iparticle + 1) = id_start(iparticle) + assigned_num(iparticle)
        end do

        !  MCMC sampling from the updated distribution
!$omp parallel do private( &
!$omp iparticle, istart, nassigned, idim, particle_cur, likelihood_cur, &
!$omp prior_cur, post_cur, jparticle, st_rand, particle_cand, &
!$omp likelihood_cand, prior_cand, post_cand, metropolis, gsvec, lsvec)
        do iparticle = 1, nparticle
            istart = id_start(iparticle)
            nassigned = assigned_num(iparticle)
            ! ----_cur means current configuration
            do idim = 1, ndim
                particle_cur(idim) = particles(idim, iparticle)
            end do
            likelihood_cur = likelihood_ls(iparticle)
            prior_cur = prior_ls(iparticle)
            !   negative log of posterior
            post_cur = gamma*likelihood_cur + prior_cur

            do jparticle = istart, istart + nassigned - 1
                ! propose particle_cand
                do idim = 1, ndim
                    st_rand(idim) = st_rand_ls(idim, jparticle)
                end do
                call dtrmv('l', 'n', 'n', ndim, cov, ndim, st_rand, 1)
                do idim = 1, ndim
                    particle_cand(idim) = particle_cur(idim) + st_rand(idim)
                    !   non negative constraints
                    ! particle_cand(idim) = max(0d0, particle_cand(idim))
                    if (particle_cand(idim) < 0d0) then
                        particle_cand(idim) = -particle_cand(idim)
                    end if
                    !   max slip constraints
                    ! particle_cand(idim) = min(max_slip, particle_cand(idim))
                    if (particle_cand(idim) > max_slip) then
                        particle_cand(idim) = 2*max_slip - particle_cand(idim)
                    end if
                end do

                ! calculate negative log likelihood/prior/posterior
                ! of the proposed configuration
                likelihood_cand = slip_calc_likelihood( &
                                  particle_cand, dvec, sigma2_full, gmat, &
                                  log_sigma_sar2, log_sigma_gnss2, nsar, ngnss, &
                                  gsvec, nobs, ndof)
                prior_cand = slip_calc_prior(particle_cand, alpha2_full, theta, nplane, nxi, neta, &
                                             lmat_index, lmat_val, lsvec, nnode)
                post_cand = gamma*likelihood_cand + prior_cand

                !  metropolis test
                metropolis = metropolis_ls(jparticle)
                if (exp(post_cur - post_cand) > metropolis) then
                    ! accept
                    do idim = 1, ndim
                        particle_cur(idim) = particle_cand(idim)
                    end do
                    likelihood_cur = likelihood_cand
                    prior_cur = prior_cand
                    post_cur = post_cand
                else
                    ! reject
                end if

                !  save to new particle list
                do idim = 1, ndim
                    particles_new(idim, jparticle) = particle_cur(idim)
                end do
                likelihood_ls_new(jparticle) = likelihood_cur
                prior_ls_new(jparticle) = prior_cur
            end do
        end do
!$omp end parallel do

        !   update configurations
!$omp parallel do private(iparticle, idim)
        do iparticle = 1, nparticle
            do idim = 1, ndim
                particles(idim, iparticle) = particles_new(idim, iparticle)
            end do
            likelihood_ls(iparticle) = likelihood_ls_new(iparticle)
            prior_ls(iparticle) = prior_ls_new(iparticle)
        end do
!$omp end parallel do
    end subroutine

    subroutine slip_smc_exec(particles, particles_new, &
                             nparticle, ndim, flag_output, &
                             output_path, llmat, &
                             max_slip, dvec, &
                             sigma2_full, alpha2_full, theta, nplane, &
                             nxi, neta, gmat, &
                             log_sigma_sar2, log_sigma_gnss2, nsar, &
                             ngnss, nobs, ndof, lmat_index, &
                             lmat_val, nnode, likelihood_ls, &
                             prior_ls, weights, mean, &
                             cov, likelihood_ls_new, &
                             prior_ls_new, assigned_num, id_start, &
                             id_dof, st_rand_ls, metropolis_ls, &
                             gsvec, lsvec, particle_cur, particle_cand, &
                             st_rand, neglog_ret)
        implicit none
        double precision, intent(in) ::  max_slip, dvec(:), &
            lmat_val(:, :), llmat(:, :), gmat(:, :), sigma2_full(:), alpha2_full(:), &
            log_sigma_sar2, log_sigma_gnss2, theta(:)
        integer, intent(in) :: nplane, nnode, ndof, nsar, ngnss, nobs, &
                               nparticle, ndim, flag_output, id_dof(:), &
                               lmat_index(:, :), nxi, neta
        character(*), intent(in) :: output_path
        double precision, intent(inout) ::  particles(:, :), &
            particles_new(:, :), likelihood_ls(:), &
            prior_ls(:), weights(:), mean(:), cov(:, :), &
            likelihood_ls_new(:), prior_ls_new(:), &
            st_rand_ls(:, :), metropolis_ls(:), gsvec(:), lsvec(:), &
            particle_cur(:), particle_cand(:), st_rand(:), neglog_ret
        integer, intent(inout) ::  assigned_num(:), id_start(:)
        integer :: iparticle, idim, idof, inode, idirection
        integer :: iter
        double precision :: gamma, neglog_evidence
        double precision, allocatable :: slip(:, :)
        iter = 0

        ! sampling from the prior distribution
        call slip_gen_init_particles(particles, likelihood_ls, prior_ls, nparticle, &
                                     ndim, llmat, max_slip, dvec, &
                                     sigma2_full, alpha2_full, theta, nplane, nxi, neta, &
                                     gmat, log_sigma_sar2, &
                                     log_sigma_gnss2, nsar, ngnss, nobs, ndof, &
                                     lmat_index, lmat_val, nnode, particle_cur, gsvec, lsvec)
        ! output result of stage 0(disabled)
        ! write (iter_char, "(i0)") iter
        ! filename = trim(trim(output_dir)//trim(iter_char)//".csv")
        ! open (17, file=filename, status='replace')
        ! do iparticle = 1, nparticle
        !    do idim = 1, ndim
        !       write (17, "(f12.5)", advance="no") particles(idim, iparticle)
        !    end do
        !    write (17, "(f12.5)") likelihood_ls(iparticle)
        ! end do
        ! close (17)
        iter = iter + 1

        gamma = 0d0
        ! product of S_j (sum for negative log value)
        neglog_ret = 0d0
        do while (1d0 - gamma > 10d-8)
            ! S_j

            ! find the gamma such that c.o.v of weights = 0.5
            gamma = slip_find_next_gamma(gamma, likelihood_ls, weights, &
                                         neglog_evidence, nparticle)
            print *, "gamma", gamma
            neglog_ret = neglog_ret + neglog_evidence

            ! normalize weights(sum of weights needs to be 1)
            call slip_normalize_weights(weights, nparticle)

            ! calculate mean and covariance of the samples
            call slip_calc_mean_particles(particles, weights, &
                                          mean, nparticle, ndim)
            call slip_calc_cov_particles(particles, weights, mean, cov, &
                                         nparticle, ndim)
            call slip_resample_particles(nparticle, weights, assigned_num)
            call slip_mcmc_sampling(gamma, particles, particles_new, &
                                    likelihood_ls, likelihood_ls_new, prior_ls, &
                                    prior_ls_new, cov, assigned_num, id_start, &
                                    nparticle, ndim, dvec, sigma2_full, alpha2_full, &
                                    theta, nplane, nxi, neta, &
                                    gmat, log_sigma_sar2, log_sigma_gnss2, nsar, &
                                    ngnss, nobs, ndof, lmat_index, lmat_val, &
                                    nnode, max_slip, st_rand_ls, metropolis_ls, &
                                    particle_cur, particle_cand, st_rand, gsvec, lsvec)
            ! output result of stage 0(disabled)
            ! write (iter_char, "(i0)") iter
            ! filename = trim(trim(output_dir)//trim(iter_char)//".csv")
            ! open (17, file=filename, status='replace')
            ! do iparticle = 1, nparticle
            !    do idim = 1, ndim
            !       write (17, "(f12.5)", advance="no") particles(idim, iparticle)
            !    end do
            !    write (17, "(f12.5)") likelihood_ls(iparticle)
            ! end do
            ! close (17)
            iter = iter + 1

        end do

        if (flag_output == 1) then
            allocate (slip(2, nnode))
            open (17, file=output_path, status='replace')
            do inode = 1, nnode
                do idirection = 1, 2
                    slip(idirection, inode) = 0d0
                end do
            end do
            do iparticle = 1, nparticle
                do idim = 1, ndim
                    particle_cur(idim) = particles(idim, iparticle)
                end do
                do idof = 1, ndof
                    inode = id_dof(idof)
                    slip(1, inode) = particle_cur(2*(idof - 1) + 1)
                    slip(2, inode) = particle_cur(2*(idof - 1) + 2)
                end do
                do inode = 1, nnode
                    do idirection = 1, 2
                        write (17, "(f12.5)", advance="no") slip(idirection, inode)
                    end do
                end do
                write (17, "(f12.5)") likelihood_ls(iparticle) + prior_ls(iparticle)
            end do
            close (17)
            deallocate (slip)
        end if

        neglog_ret = min(1d10, neglog_ret)

    end subroutine
end module smc_slip
