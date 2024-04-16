module smc_slip
    use, intrinsic :: ieee_arithmetic
    use type_mat
    use omp_lib
    implicit none
contains
    subroutine random_number_correction(x)
        implicit none
        double precision, intent(inout) :: x
        call random_number(x)
        if (x < 1d-5) then
            x = 1d-5
        end if
    end subroutine random_number_correction

    subroutine slip_BoxMuller(p, q)
        implicit none
        double precision, intent(inout) :: p, q
        double precision :: pi, r, s
        pi = 2d0*asin(1d0)
        !   uniform random numbers between 0 and 1
        call random_number_correction(r)
        call random_number_correction(s)
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
                                                   ngnss, gsvec, &
                                                   nobs, ndof, gsgmat, &
                                                   gsdvec, dsd, nsar_total)
        implicit none
        double precision, intent(in) :: svec(:), dvec(:), sigma2_full(:), &
            gmat(:, :), log_sigma_sar2, log_sigma_gnss2, gsgmat(:, :), &
            gsdvec(:), dsd
        integer, intent(in) :: ngnss, nobs, ndof, nsar_total
        double precision, intent(inout) :: gsvec(:)
        integer :: i, j

        double precision :: dmax, st_time, en_time

        slip_calc_likelihood = 0d0
        do i = 1, 2*ndof
            do j = 1, i - 1
                slip_calc_likelihood = slip_calc_likelihood &
                                        + svec(i) * gsgmat(j, i) * svec(j)
            end do
        end do
        slip_calc_likelihood = slip_calc_likelihood * 2d0
        do i = 1, 2*ndof
            slip_calc_likelihood = slip_calc_likelihood &
                                   + svec(i) * gsgmat(i, i) * svec(i)
        end do
        slip_calc_likelihood = (slip_calc_likelihood &
                                - 2d0*dot_product(gsdvec, svec) &
                                + dsd)/2d0 &
                                + nsar_total*log_sigma_sar2/2d0 + 3d0*ngnss*log_sigma_gnss2/2d0

        ! do i = 1, nobs
        !     gsvec(i) = 0d0
        ! end do
        ! call dgemv('n', nobs, 2*ndof, 1d0, gmat, &
        !            nobs, svec, 1, 0d0, gsvec, 1)

        ! slip_calc_likelihood = nsar_total*log_sigma_sar2/2d0 + 3d0*ngnss*log_sigma_gnss2/2d0
        ! ! slip_calc_likelihood = 0d0
        ! do i = 1, nobs
        !     ! slip_calc_likelihood = slip_calc_likelihood + & 
        !     !                         log(sigma2_full(i)) / 2d0
        !     slip_calc_likelihood = slip_calc_likelihood + &
        !                            (dvec(i) - gsvec(i))**2d0/(2d0*sigma2_full(i))
        ! end do
        ! print *, slip_calc_likelihood
        ! stop
    end function slip_calc_likelihood

    double precision function slip_calc_prior(svec, alpha2_full, theta, nplane, &
                                              lmat_index, lmat_val, lsvec, nnode)
        implicit none
        double precision, intent(in) :: svec(:), alpha2_full(:), &
            lmat_val(:, :), theta(:)
        integer, intent(in) :: lmat_index(:, :), nnode, nplane
        double precision, intent(inout) :: lsvec(:)
        integer :: n, i, j
        double precision :: tmp, log_alpha2

        n = 2*nnode
        do i = 1, n
            lsvec(i) = 0d0
        end do
        do i = 1, n
            tmp = 0d0
            do j = 1, 9
                tmp = tmp + lmat_val(j, i)*svec(lmat_index(j, i))
            end do
            lsvec(i) = tmp
        end do

        slip_calc_prior = 0d0
        ! do i = 1, nplane
        !     log_alpha2 = theta(8*i)
        !     slip_calc_prior = slip_calc_prior + &
        !                       2d0*(nxi + 1)*(neta + 1)*log_alpha2/2d0
        ! end do
        do i = 1, n
            slip_calc_prior = slip_calc_prior + lsvec(i)**2d0/(2d0*alpha2_full(i))
        end do

    end function slip_calc_prior

    subroutine slip_calc_grad(grad, svec, gmat, lmat_index, lmat_val, &
                              ltmat_index, ltmat_val, dvec, gamma, &
                              gsvec, lsvec, sigma2_full, alpha2_full, &
                              nobs, ndof, nnode, ndim, gsdvec, gsgmat)
        implicit none
        double precision, intent(in) :: svec(:), gmat(:, :), dvec(:), gamma, &
            sigma2_full(:), alpha2_full(:), lmat_val(:, :), ltmat_val(:, :), &
            gsdvec(:), gsgmat(:, :)
        integer, intent(in) :: nobs, ndof, nnode, ndim, &
                               lmat_index(:, :), ltmat_index(:, :)
        double precision, intent(inout) :: grad(:), gsvec(:), lsvec(:)
        double precision :: grad1(ndim), grad2(ndim), tmp, tmp_vec(nobs)
        integer :: inode, iobs, idof, idim, j

        grad = 0d0

        ! grad1 = 0d0
        ! call dgemv('n', nobs, 2*ndof, 1d0, gmat, &
        !            nobs, svec, 1, 0d0, gsvec, 1)
        ! do iobs = 1, nobs
        !     tmp_vec(iobs) = (dvec(iobs) - gsvec(iobs))/sigma2_full(iobs)
        ! end do
        ! call dgemv('t', nobs, 2*ndof, -gamma, gmat, &
        !            nobs, tmp_vec, 1, 0d0, grad1, 1)

        grad1 = 0d0
        call dgemv('n', ndim, ndim, 1d0, gsgmat, ndim, &
                   svec, 1, 0d0, grad1, 1)
        do idim = 1, ndim
            grad1(idim) = grad1(idim) - gsdvec(idim)
            grad1(idim) = grad1(idim)*gamma
        end do

        do inode = 1, 2*nnode
            tmp = 0d0
            do j = 1, 9
                tmp = tmp + lmat_val(j, inode)*svec(lmat_index(j, inode))
            end do
            lsvec(inode) = tmp/alpha2_full(inode)
        end do

        do idof = 1, 2*ndof
            tmp = 0d0
            do j = 1, 9
                tmp = tmp + ltmat_val(j, idof)*lsvec(ltmat_index(j, idof))
            end do
            grad2(idof) = tmp
        end do

        do idim = 1, ndim
            grad(idim) = grad1(idim) + grad2(idim)
        end do
    end subroutine slip_calc_grad

    subroutine slip_gen_init_particles(particles, likelihood_ls, prior_ls, nparticle, &
                                       ndim, llmat, max_slip, dvec, &
                                       sigma2_full, alpha2_full, theta, nplane, &
                                       gmat, log_sigma_sar2, &
                                       log_sigma_gnss2, ngnss, nobs, ndof, &
                                       lmat_index, lmat_val, nnode, particle_cur, gsvec, lsvec, &
                                       gsdvec, dsd, nsar_total, gsgmat)
        implicit none
        double precision, intent(in) :: max_slip, dvec(:), &
            lmat_val(:, :), llmat(:, :), gmat(:, :), sigma2_full(:), &
            alpha2_full(:), log_sigma_sar2, log_sigma_gnss2, theta(:), &
            gsdvec(:), dsd, gsgmat(:,:)
        integer, intent(in) :: nnode, ndof, ngnss, nobs, &
                               nparticle, ndim, lmat_index(:, :), nplane, &
                               nsar_total
        double precision, intent(inout) ::  particles(:, :), &
            likelihood_ls(:), prior_ls(:), particle_cur(:), gsvec(:), lsvec(:)
        integer :: iparticle, idim, jdim, cnt
        double precision :: x, fx, f1, f0
        double precision :: y_i, y_prev, err
        double precision :: v1, v2  ! box muller
        double precision :: sigma2_i, mu_i
        double precision :: st_time, en_time

        st_time = omp_get_wtime()
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
                if (mod(idim, 2) == 0) then
                    ! dip slip
                    f0 = cdf_norm(0d0, mu_i, sigma2_i)
                else
                    ! strike slip
                    f0 = cdf_norm(-max_slip, mu_i, sigma2_i)
                end if
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
                ! y_i = min(y_i, max_slip)
                ! y_i = max(y_i, 0d0)
                particle_cur(idim) = y_i
            end do
            do idim = 1, ndim
                particles(idim, iparticle) = particle_cur(idim)
            end do
        end do
        en_time = omp_get_wtime()

!$omp parallel do private(iparticle, idim, particle_cur, gsvec, lsvec)
        do iparticle = 1, nparticle
            do idim = 1, ndim
                particle_cur(idim) = particles(idim, iparticle)
            end do
            !   calculate negative log likelihood and prior
            likelihood_ls(iparticle) = &
                slip_calc_likelihood(particle_cur, dvec, sigma2_full, &
                                     gmat, log_sigma_sar2, log_sigma_gnss2, &
                                     ngnss, gsvec, nobs, ndof, gsgmat, gsdvec, dsd, nsar_total)
            prior_ls(iparticle) = slip_calc_prior(particle_cur, alpha2_full, theta, nplane, &
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
        integer :: iparticle, cnt
        double precision :: min_likelihood, cv_threshold, lower, upper, &
            err, gamma, diff_gamma, mean, std, cv, likelihood, evidence

        ! find minimum of negative log likelihood
        min_likelihood = likelihood_ls(1)
!$omp parallel do private(iparticle) reduction(min : min_likelihood)
        do iparticle = 1, nparticle
            min_likelihood = min(min_likelihood, likelihood_ls(iparticle))
        end do
!$omp end parallel do

        cv_threshold = 2d0
        ! binary search for the next gamma
        ! such that c.o.v of the weight is equivalent to cv_threashold
        lower = gamma_prev
        upper = 1d0
        err = 1d0
        gamma = 1d0
        cnt = 0
        do while (err > 1d-10)
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
            if (abs(gamma - 1d0) < 1d-10) then
                exit
            end if
            ! if (cnt > 1000) then
            !     exit
            ! end if
            cnt = cnt + 1
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

    double precision function slip_find_next_gamma_ess(gamma_prev, likelihood_ls, &
                                                       weights, neglog_evidence, &
                                                       nparticle)
        implicit none
        double precision, intent(in) :: gamma_prev, likelihood_ls(:)
        integer, intent(in) :: nparticle
        double precision, intent(inout) :: neglog_evidence, weights(:)
        integer :: iparticle
        double precision :: min_likelihood, cv_threshold, lower, upper, &
            err, gamma, diff_gamma, mean, std, cv, likelihood, evidence
        double precision :: wsum, w2sum, ess, ess_threshold

        ess_threshold = nparticle/2d0
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

            ! calculate effective sample size
            wsum = 0d0; 
            w2sum = 0d0; 
!$omp parallel do private(iparticle) reduction(+:wsum, w2sum)
            do iparticle = 1, nparticle
                wsum = wsum + weights(iparticle)
                w2sum = w2sum + weights(iparticle)**2
            end do
            ess = wsum**2/w2sum

            if (ess < ess_threshold) then
                upper = gamma
            else
                lower = gamma
            end if
            err = abs(ess - ess_threshold)
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
        slip_find_next_gamma_ess = gamma
    end function slip_find_next_gamma_ess

    subroutine slip_normalize_weights(weights, nparticle)
        implicit none
        double precision, intent(inout) :: weights(:)
        integer, intent(in) :: nparticle
        double precision :: sum
        integer :: iparticle

        sum = 0d0
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
                                       cov, nparticle, ndim, cov_diag)
        implicit none
        double precision, intent(in) :: particles(:, :), weights(:), mean(:)
        integer, intent(in) :: nparticle, ndim
        double precision, intent(inout) :: cov(:, :), cov_diag(:)
        integer :: iparticle, idim, jdim, ierr
        double precision weight, di, dj
!         do jdim = 1, ndim
!             do idim = jdim, ndim
!                 cov(idim, jdim) = 0d0
!             end do
!         end do

! !$omp parallel do private(iparticle, weight, idim, di, jdim, dj) &
! !$omp reduction(+:cov)
!         do iparticle = 1, nparticle
!             weight = weights(iparticle)
!             do jdim = 1, ndim
!                 dj = (particles(jdim, iparticle) - mean(jdim))
!                 do idim = jdim, ndim
!                     di = (particles(idim, iparticle) - mean(idim))
!                     cov(idim, jdim) = cov(idim, jdim) + weight*di*dj
!                 end do
!             end do
!         end do
! !$omp end parallel do
!         do idim = 1, ndim
!             cov_diag(idim) = cov(idim, idim)
!         end do
!         do jdim = 1, ndim
!             do idim = jdim, ndim
!                 cov(idim, jdim) = cov(idim, jdim)*0.04d0
!             end do
!         end do

!         ! LAPACK function for LU decomposition of matrix
!         call dpotrf('L', ndim, cov, ndim, ierr)

        cov_diag = 0d0
!$omp parallel do private(iparticle, weight, idim, di) &
!$omp reduction(+:cov_diag)
        do iparticle = 1, nparticle
            weight = weights(iparticle)
            do idim = 1, ndim
                di = particles(idim, iparticle) - mean(idim)
                cov_diag(idim) = cov_diag(idim) + weight*di*di
            end do
        end do
!$omp end parallel do

        do idim = 1, ndim
            if (abs(cov_diag(idim)) < 1d-3) then
                cov_diag(idim) = 1d-3
            end if
        end do

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
        call random_number_correction(u)
        u = u/d_nparticle
        do iparticle = 1, nparticle
            assigned_num(iparticle) = &
                floor((weights(iparticle) - u)*d_nparticle) + 1
            u = u + assigned_num(iparticle)/d_nparticle - weights(iparticle)
        end do
    end subroutine slip_resample_particles

    subroutine slip_hmc_sampling(gamma, particles, particles_new, &
                                 likelihood_ls, likelihood_ls_new, prior_ls, &
                                 prior_ls_new, cov, cov_diag, assigned_num, id_start, &
                                 nparticle, ndim, dvec, sigma2_full, alpha2_full, &
                                 theta, nplane, &
                                 gmat, log_sigma_sar2, log_sigma_gnss2, &
                                 ngnss, nobs, ndof, lmat_index, lmat_val, &
                                 ltmat_index, ltmat_val, &
                                 nnode, max_slip, st_rand_ls, metropolis_ls, &
                                 particle_cur, particle_cand, st_rand, gsvec, lsvec, &
                                 dtau_ls, ntau_ls, gsdvec, gsgmat, ntau_upper,  &
                                 dsd, nsar_total)
        implicit none
        double precision, intent(in) :: max_slip, dvec(:), &
            lmat_val(:, :), ltmat_val(:, :), &
            gmat(:, :), sigma2_full(:), alpha2_full(:), theta(:), &
            log_sigma_sar2, log_sigma_gnss2, gamma, cov(:, :), cov_diag(:), &
            dtau_ls(:), gsdvec(:), gsgmat(:, :), dsd
        integer, intent(in) :: nnode, ndof, ngnss, nobs, nplane,  &
                               nparticle, ndim, lmat_index(:, :), ltmat_index(:, :), &
                               assigned_num(:), ntau_ls(:), ntau_upper, nsar_total
        double precision, intent(inout) ::  particles(:, :), &
            particles_new(:, :), likelihood_ls(:), prior_ls(:), &
            likelihood_ls_new(:), prior_ls_new(:), &
            st_rand_ls(:, :), metropolis_ls(:), particle_cur(:), &
            particle_cand(:), st_rand(:), gsvec(:), lsvec(:)

        integer, intent(inout) :: id_start(:)
        integer :: iparticle, jparticle, idim, nassigned, istart, ierr
        double precision ::likelihood_cur, likelihood_cand, &
            prior_cur, prior_cand, post_cur, post_cand, metropolis
        double precision :: v1, v2 ! box muller

        ! for HMC
        double precision :: dtau, pvec(ndim), grad(2*ndof), ham_cur, ham_cand
        double precision :: acc_rate, gamma_tmp, delta
        integer :: itau, ntau

        ! timer
        double precision :: st_time, en_time

        st_time = omp_get_wtime()
        do iparticle = 1, nparticle
            do idim = 1, ndim
                call slip_BoxMuller(v1, v2)
                st_rand_ls(idim, iparticle) = v1
            end do
            call random_number_correction(metropolis)
            metropolis_ls(iparticle) = metropolis
        end do
        en_time = omp_get_wtime()

        id_start(1) = 1
        do iparticle = 1, nparticle - 1
            id_start(iparticle + 1) = id_start(iparticle) + assigned_num(iparticle)
        end do

        acc_rate = 0d0
!$omp parallel do private(&
!$omp iparticle, istart, nassigned, idim, particle_cur, likelihood_cur, &
!$omp prior_cur, post_cur, jparticle, particle_cand, st_rand, pvec, &
!$omp ham_cur, itau, grad, gsvec, lsvec, likelihood_cand, prior_cand, &
!$omp post_cand, ham_cand, metropolis, dtau, ntau) reduction(+:acc_rate)
        ! hamiltonian monte carlo
        do iparticle = 1, nparticle
            istart = id_start(iparticle)
            nassigned = assigned_num(iparticle)
            ! ----_cur means current configuration
            do idim = 1, ndim
                particle_cur(idim) = particles(idim, iparticle)
            end do
            likelihood_cur = likelihood_ls(iparticle)
            prior_cur = prior_ls(iparticle)
            ! negative log of posterior
            post_cur = gamma*likelihood_cur + prior_cur

            do jparticle = istart, istart + nassigned - 1
                dtau = dtau_ls(jparticle)
                ntau = ntau_ls(jparticle)
                ntau = 50
                ! if (ntau > ntau_upper) then
                !     print *, "ntau: ", ntau, " is too large"
                !     print *, "theta: ", particle_cur
                !     stop
                ! end if
                particle_cand = particle_cur
                ! sampling momentum
                do idim = 1, ndim
                    st_rand(idim) = st_rand_ls(idim, jparticle)
                    pvec(idim) = st_rand(idim)/sqrt(cov_diag(idim))
                    ! if (ieee_is_nan(pvec(idim))) then
                    !     print *, "pvec is nan"
                    !     print *, "st_rand: ",st_rand
                    !     print *, "cov_diag: ", cov_diag
                    !     stop
                    ! end if
                end do
                call slip_leapfrog(ham_cur, post_cur, ndim, pvec, cov_diag, &
                                   particle_cand, dtau, ntau, grad, gmat, lmat_index, lmat_val, &
                                   ltmat_index, ltmat_val, dvec, gamma, gsvec, lsvec, &
                                   sigma2_full, alpha2_full, nobs, ndof, nnode, max_slip, &
                                   likelihood_cand, log_sigma_sar2, log_sigma_gnss2, &
                                   ngnss, prior_cand, theta, nplane, &
                                   post_cand, ham_cand, particle_cur, gsdvec, gsgmat, &
                                   dsd, nsar_total)
                !  metropolis test
                metropolis = metropolis_ls(jparticle)
                if (((ham_cur - ham_cand) > log(metropolis))) then
                    ! accept
                    do idim = 1, ndim
                        particle_cur(idim) = particle_cand(idim)
                    end do
                    likelihood_cur = likelihood_cand
                    prior_cur = prior_cand
                    post_cur = post_cand
                    acc_rate = acc_rate + 1d0/nparticle
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
        ! print *, "acc_rate: ", acc_rate
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
    end subroutine slip_hmc_sampling

    subroutine slip_leapfrog(ham_cur, post_cur, ndim, pvec, cov_diag, &
                             particle_cand, dtau, ntau, grad, gmat, lmat_index, lmat_val, &
                             ltmat_index, ltmat_val, dvec, gamma, gsvec, lsvec, &
                             sigma2_full, alpha2_full, nobs, ndof, nnode, max_slip, &
                             likelihood_cand, log_sigma_sar2, log_sigma_gnss2, &
                             ngnss, prior_cand, theta, nplane, &
                             post_cand, ham_cand, particle_cur, gsdvec, gsgmat, &
                             dsd, nsar_total)
        implicit none
        integer, intent(in) :: ndim, ntau, lmat_index(:, :), ltmat_index(:, :), &
                               nobs, ndof, nnode, ngnss, nplane,  &
                               nsar_total
        double precision, intent(in) :: post_cur, cov_diag(:), dtau, &
            gmat(:, :), lmat_val(:, :), ltmat_val(:, :), dvec(:), gamma, &
            sigma2_full(:), alpha2_full(:), max_slip, log_sigma_sar2, &
            log_sigma_gnss2, theta(:), particle_cur(:), gsdvec(:), gsgmat(:, :), dsd
        double precision, intent(inout) :: ham_cur, pvec(:), particle_cand(:), grad(:), &
            gsvec(:), lsvec(:), likelihood_cand, prior_cand, post_cand, ham_cand
        integer :: idim, itau

        double precision :: dtmp1, dtmp2

        ! initial hamiltonian
        ham_cur = post_cur
        do idim = 1, ndim
            ham_cur = ham_cur + pvec(idim)*cov_diag(idim)*pvec(idim)/2d0
        end do
        ! leapfrog integration
        do idim = 1, ndim
            particle_cand(idim) = particle_cand(idim) &
                                  + cov_diag(idim)*pvec(idim)*5d-1*dtau
            if (abs(particle_cand(idim)) > 1d3) then
                particle_cand = particle_cur
                ham_cand = 1d20
                return
            end if
            ! reflection
            if (mod(idim, 2) == 0) then
                ! dip slip
                do while (particle_cand(idim) < 0d0 .or. &
                        particle_cand(idim) > max_slip)
                    if (particle_cand(idim) > max_slip) then
                        particle_cand(idim) = &
                            2d0*max_slip - particle_cand(idim)
                        pvec(idim) = -pvec(idim)
                    end if
                    if (particle_cand(idim) < 0d0) then
                        particle_cand(idim) = &
                            -particle_cand(idim)
                        pvec(idim) = -pvec(idim)
                    end if
                end do
            else
                ! strike slip
                do while (particle_cand(idim) < -max_slip .or. &
                        particle_cand(idim) > max_slip)
                    if (particle_cand(idim) > max_slip) then
                        particle_cand(idim) = &
                            2d0*max_slip - particle_cand(idim)
                        pvec(idim) = -pvec(idim)
                    end if
                    if (particle_cand(idim) < -max_slip) then
                        particle_cand(idim) = &
                            -2d0*max_slip -particle_cand(idim)
                        pvec(idim) = -pvec(idim)
                    end if
                end do
            end if
            ! if (ieee_is_nan(particle_cand(idim))) then
            !     print *, "particle_cand is nan"
            !     print *, "pvec: ", pvec
            !     print *, "cov_diag: ", cov_diag
            !     print *, "dtau: ", dtau
            !     stop
            ! end if
        end do
        do itau = 1, ntau
            call slip_calc_grad(grad, particle_cand, gmat, lmat_index, lmat_val, &
                                ltmat_index, ltmat_val, dvec, 5d-1, &
                                gsvec, lsvec, sigma2_full, alpha2_full, &
                                nobs, ndof, nnode, ndim, gsdvec, gsgmat)
            
            ! print *, "check gradient"
            ! open(10, file="tmp/svec", status="replace")
            ! write(10, *), particle_cand
            ! close(10)
            ! print *, "svec: ", particle_cand
            ! print *, "-----------------"
            ! do idim = 1, ndim
            !     particle_cand(idim) = particle_cand(idim) - 1d-8
            !     dtmp1 = 5d-1*slip_calc_likelihood(particle_cand, dvec, sigma2_full, gmat, log_sigma_sar2, log_sigma_gnss2, ngnss, gsvec, nobs, ndof, gsgmat_l, gsdvec, dsd, nsar_total)
            !     dtmp1 = dtmp1 + slip_calc_prior(particle_cand, alpha2_full, theta, nplane, lmat_index, lmat_val, lsvec, nnode)
            !     particle_cand(idim) = particle_cand(idim) + 2d-8
            !     dtmp2 = 5d-1*slip_calc_likelihood(particle_cand, dvec, sigma2_full, gmat, log_sigma_sar2, log_sigma_gnss2, ngnss, gsvec, nobs, ndof, gsgmat_l, gsdvec, dsd, nsar_total)
            !     dtmp2 = dtmp2 + slip_calc_prior(particle_cand, alpha2_full, theta, nplane, lmat_index, lmat_val, lsvec, nnode)
            !     particle_cand(idim) = particle_cand(idim) - 1d-8
            !     print *, (dtmp2 - dtmp1)/2d-8, grad(idim)
            ! end do
            ! stop
            call slip_calc_grad(grad, particle_cand, gmat, lmat_index, lmat_val, &
                                ltmat_index, ltmat_val, dvec, gamma, &
                                gsvec, lsvec, sigma2_full, alpha2_full, &
                                nobs, ndof, nnode, ndim, gsdvec, gsgmat)
            do idim = 1, ndim
                pvec(idim) = pvec(idim) - grad(idim)*dtau
                ! if (ieee_is_nan(pvec(idim))) then
                !     print *, "pvec is nan"
                !     print *, "grad: ", grad
                !     print *, "dtau: ", dtau
                !     stop
                ! end if
            end do
            do idim = 1, ndim
                particle_cand(idim) = particle_cand(idim) &
                                      + cov_diag(idim)*pvec(idim)*dtau
                if (abs(particle_cand(idim)) > 1d3) then
                    particle_cand = particle_cur
                    ham_cand = 1d20
                    return
                end if
                ! reflection
                if (mod(idim, 2) == 0) then
                    ! dip slip
                    do while (particle_cand(idim) < 0d0 .or. &
                            particle_cand(idim) > max_slip)
                        if (particle_cand(idim) > max_slip) then
                            particle_cand(idim) = &
                                2d0*max_slip - particle_cand(idim)
                            pvec(idim) = -pvec(idim)
                        end if
                        if (particle_cand(idim) < 0d0) then
                            particle_cand(idim) = &
                                -particle_cand(idim)
                            pvec(idim) = -pvec(idim)
                        end if
                    end do
                else
                    ! strike slip
                    do while (particle_cand(idim) < -max_slip .or. &
                            particle_cand(idim) > max_slip)
                        if (particle_cand(idim) > max_slip) then
                            particle_cand(idim) = &
                                2d0*max_slip - particle_cand(idim)
                            pvec(idim) = -pvec(idim)
                        end if
                        if (particle_cand(idim) < -max_slip) then
                            particle_cand(idim) = &
                                -2d0*max_slip -particle_cand(idim)
                            pvec(idim) = -pvec(idim)
                        end if
                    end do
                end if
                ! if (ieee_is_nan(particle_cand(idim))) then
                !     print *, "particle_cand is nan"
                !     print *, "pvec: ", pvec
                !     print *, "cov_diag: ", cov_diag
                !     print *, "dtau: ", dtau
                !     stop
                ! end if
            end do
        end do
        call slip_calc_grad(grad, particle_cand, gmat, lmat_index, lmat_val, &
                            ltmat_index, ltmat_val, dvec, gamma, &
                            gsvec, lsvec, sigma2_full, alpha2_full, &
                            nobs, ndof, nnode, ndim, gsdvec, gsgmat)
        do idim = 1, ndim
            pvec(idim) = pvec(idim) - grad(idim)*5d-1*dtau
            ! if (ieee_is_nan(pvec(idim))) then
            !     print *, "pvec is nan"
            !     print *, "grad: ", grad
            !     print *, "dtau: ", dtau
            !     stop
            ! end if
        end do
        ! do idim = 1, ndim
        !     !   non negative constraints
        !     if (particle_cand(idim) < 0d0) then
        !         particle_cand(idim) = min(-particle_cand(idim), max_slip)
        !         ! particle_cand = particle_cur
        !         ! ham_cand = 1d20
        !         ! return
        !     end if
        !     !   max slip constraints
        !     if (particle_cand(idim) > max_slip) then
        !         particle_cand(idim) = max(2*max_slip - particle_cand(idim), 0d0)
        !         ! particle_cand = particle_cur
        !         ! ham_cand = 1d20
        !         ! return
        !     end if
        ! end do
        ! final hamiltonian
        likelihood_cand = slip_calc_likelihood( &
                          particle_cand, dvec, sigma2_full, gmat, &
                          log_sigma_sar2, log_sigma_gnss2, ngnss, &
                          gsvec, nobs, ndof, gsgmat, gsdvec, dsd, nsar_total)
        prior_cand = slip_calc_prior(particle_cand, alpha2_full, theta, nplane, &
                                     lmat_index, lmat_val, lsvec, nnode)
        post_cand = gamma*likelihood_cand + prior_cand
        ham_cand = post_cand
        do idim = 1, ndim
            ham_cand = ham_cand + pvec(idim)*cov_diag(idim)*pvec(idim)/2d0
        end do
    end subroutine

    subroutine slip_hmc_tuning(gamma, particles, particles_new, &
                               likelihood_ls, likelihood_ls_new, prior_ls, &
                               prior_ls_new, cov, cov_diag, assigned_num, id_start, &
                               nparticle, ndim, dvec, sigma2_full, alpha2_full, &
                               theta, nplane, &
                               gmat, log_sigma_sar2, log_sigma_gnss2, &
                               ngnss, nobs, ndof, lmat_index, lmat_val, &
                               ltmat_index, ltmat_val, &
                               nnode, max_slip, st_rand_ls, metropolis_ls, &
                               particle_cur, particle_cand, st_rand, gsvec, lsvec, &
                               log_dtau_upper, log_dtau_lower, ntau_upper, ntau_lower, &
                               dtau_ls, ntau_ls, &
                               gsdvec, gsgmat, tuning_factor, &
                               dsd, nsar_total)
        implicit none
        double precision, intent(in) :: max_slip, dvec(:), &
            lmat_val(:, :), ltmat_val(:, :), &
            gmat(:, :), sigma2_full(:), alpha2_full(:), theta(:), &
            log_sigma_sar2, log_sigma_gnss2, gamma, cov(:, :), cov_diag(:), &
            gsdvec(:), gsgmat(:, :), dsd
        integer, intent(in) :: nnode, ndof, ngnss, nobs, nplane,  &
                               nparticle, ndim, lmat_index(:, :), ltmat_index(:, :), &
                               assigned_num(:), tuning_factor, nsar_total
        double precision, intent(inout) ::  particles(:, :), &
            particles_new(:, :), likelihood_ls(:), prior_ls(:), &
            likelihood_ls_new(:), prior_ls_new(:), &
            st_rand_ls(:, :), metropolis_ls(:), particle_cur(:), &
            particle_cand(:), st_rand(:), gsvec(:), lsvec(:), dtau_ls(:), &
            log_dtau_upper, log_dtau_lower
        

        integer, intent(inout) :: id_start(:), ntau_ls(:), ntau_upper, ntau_lower
        integer :: iparticle, jparticle, kparticle, idim, nassigned, istart, ierr
        double precision ::likelihood_cur, likelihood_cand, &
            prior_cur, prior_cand, post_cur, post_cand, metropolis
        double precision :: v1, v2 ! box muller

        ! for HMC
        double precision :: log_dtau, dtau, pvec(ndim), grad(2*ndof), ham_cur, ham_cand
        double precision :: acc_rate, gamma_tmp, delta
        integer :: itau, ntau

        ! for tuning HMC
        double precision :: weights_hmc(nparticle), sum_score, dtau_mean
        double precision :: dham, score, d_nparticle, u, &
            dham_ls(nparticle), score_ls(nparticle)
        double precision :: dtau_ls_new(nparticle)
        integer :: iassigned, ntau_ls_new(nparticle), &
                   nassigned_hmc(nparticle), cnt, ntau_mean

        ! timer
        double precision :: st_time, en_time

        d_nparticle = dble(nparticle)
        id_start(1) = 1
        do iparticle = 1, nparticle - 1
            id_start(iparticle + 1) = id_start(iparticle) + assigned_num(iparticle)
        end do

        cnt = 0
        ! sampling from uniform distribution
        do iparticle = 1, nparticle, tuning_factor
            istart = id_start(iparticle)
            nassigned = assigned_num(iparticle)
            do jparticle = istart, istart + nassigned - 1
                cnt = cnt + 1
                call random_number_correction(dtau)
                ntau = int(ntau_lower + dtau*(ntau_upper - ntau_lower))
                ntau_ls(jparticle) = ntau
                call random_number_correction(log_dtau)
                log_dtau = log_dtau_lower + log_dtau*(log_dtau_upper - log_dtau_lower)
                dtau = exp(log_dtau)
                dtau_ls(jparticle) = dtau
            end do
        end do
        if (cnt == 0) then
            return
        end if

        st_time = omp_get_wtime()
        ! do iparticle = 1, nparticle, tuning_factor
        do iparticle = 1, nparticle
            istart = id_start(iparticle)
            nassigned = assigned_num(iparticle)
            do jparticle = istart, istart + nassigned - 1
                do idim = 1, ndim
                    call slip_BoxMuller(v1, v2)
                    st_rand_ls(idim, jparticle) = v1
                end do
                call random_number_correction(metropolis)
                metropolis_ls(jparticle) = metropolis
            end do
        end do
        en_time = omp_get_wtime()

        score_ls = -1d0
!$omp parallel do private(&
!$omp iparticle, istart, nassigned, idim, particle_cur, likelihood_cur, &
!$omp prior_cur, post_cur, jparticle, particle_cand, st_rand, pvec, &
!$omp ham_cur, itau, grad, gsvec, lsvec, likelihood_cand, prior_cand, &
!$omp post_cand, ham_cand, metropolis, dtau, ntau, dham, score) reduction(+:acc_rate)
        ! hamiltonian monte carlo
        do iparticle = 1, nparticle, tuning_factor
            istart = id_start(iparticle)
            nassigned = assigned_num(iparticle)
            ! ----_cur means current configuration
            do idim = 1, ndim
                particle_cur(idim) = particles(idim, iparticle)
            end do
            likelihood_cur = likelihood_ls(iparticle)
            prior_cur = prior_ls(iparticle)
            ! negative log of posterior
            post_cur = gamma*likelihood_cur + prior_cur

            do jparticle = istart, istart + nassigned - 1
                dtau = dtau_ls(jparticle)
                ! ntau = ntau_ls(jparticle)
                ntau = 50
                particle_cand = particle_cur
                ! sampling momentum
                do idim = 1, ndim
                    st_rand(idim) = st_rand_ls(idim, jparticle)
                    pvec(idim) = st_rand(idim)/sqrt(cov_diag(idim))
                    ! if (ieee_is_nan(pvec(idim))) then
                    !     print *, "pvec is nan"
                    !     print *, "st_rand: ",st_rand
                    !     print *, "cov_diag: ", cov_diag
                    !     stop
                    ! end if
                end do
                call slip_leapfrog(ham_cur, post_cur, ndim, pvec, cov_diag, &
                                   particle_cand, dtau, ntau, grad, gmat, lmat_index, lmat_val, &
                                   ltmat_index, ltmat_val, dvec, gamma, gsvec, lsvec, &
                                   sigma2_full, alpha2_full, nobs, ndof, nnode, max_slip, &
                                   likelihood_cand, log_sigma_sar2, log_sigma_gnss2, &
                                   ngnss, prior_cand, theta, nplane, &
                                   post_cand, ham_cand, particle_cur, gsdvec, gsgmat, &
                                   dsd, nsar_total)
                dham = (ham_cur - ham_cand)

                score = 0d0
                do idim = 1, ndim
                    score = score + ((particle_cand(idim) - particle_cur(idim))**2)/cov_diag(idim)
                end do
                ! score = score/ntau*exp(min(0d0, dham))
                score = score*exp(min(0d0, dham))
                if (ieee_is_nan(score)) then
                    print *, "score is nan"
                    print *, "particle_cand: ", particle_cand
                    print *, "particle_cur: ", particle_cur
                    print *, "pvec:", pvec
                    print *, "cov_diag: ", cov_diag
                    print *, "ham_cand: ", ham_cand
                    print *, "ham_cur: ", ham_cur
                    print *, "likelihood_cand: ", likelihood_cand
                    print *, "likelihood_cur: ", likelihood_cur
                    print *, "prior_cand: ", prior_cand
                    print *, "prior_cur: ", prior_cur
                    stop
                end if
                score_ls(jparticle) = score

                !  metropolis test
                metropolis = metropolis_ls(jparticle)
                if (ham_cur - ham_cand > log(metropolis)) then
                    ! accept
                    do idim = 1, ndim
                        particle_cur(idim) = particle_cand(idim)
                    end do
                    likelihood_cur = likelihood_cand
                    prior_cur = prior_cand
                    post_cur = post_cand
                    acc_rate = acc_rate + 1d0/dble(nparticle/tuning_factor)
                else
                    ! reject
                end if
            end do
        end do

        sum_score = 0d0
!$omp parallel do private(iparticle, istart, nassigned, jparticle) reduction(+:sum_score)
        do iparticle = 1, nparticle, tuning_factor
            istart = id_start(iparticle)
            nassigned = assigned_num(iparticle)
            do jparticle = istart, istart + nassigned - 1
                sum_score = sum_score + score_ls(jparticle)
            end do
        end do
!$omp end parallel do

        ! if degenerated
        if (sum_score < 1d-5) then
            sum_score = 0d0
!$omp parallel do private(iparticle, istart, nassigned, jparticle) reduction(+:sum_score)
            do iparticle = 1, nparticle, tuning_factor
                istart = id_start(iparticle)
                nassigned = assigned_num(iparticle)
                do jparticle = istart, istart + nassigned - 1
                    score_ls(jparticle) = 1d0
                    sum_score = sum_score + 1d0
                end do
            end do
!$omp end parallel do
        end if
        weights_hmc = -1d0
!$omp parallel do private(iparticle, istart, nassigned, jparticle)
        do iparticle = 1, nparticle, tuning_factor
            istart = id_start(iparticle)
            nassigned = assigned_num(iparticle)
            do jparticle = istart, istart + nassigned - 1
                weights_hmc(jparticle) = score_ls(jparticle)/sum_score
            end do
        end do
!$omp end parallel do

        ! residual systematic resampling
        call random_number_correction(u)
        nassigned_hmc = 0
        u = u/d_nparticle
        do iparticle = 1, nparticle, tuning_factor
            istart = id_start(iparticle)
            nassigned = assigned_num(iparticle)
            do jparticle = istart, istart + nassigned - 1
                nassigned_hmc(jparticle) = &
                    floor((weights_hmc(jparticle) - u)*d_nparticle) + 1
                u = u + nassigned_hmc(jparticle)/d_nparticle - weights_hmc(jparticle)
            end do
        end do

        cnt = 0
        do iparticle = 1, nparticle, tuning_factor
            istart = id_start(iparticle)
            nassigned = assigned_num(iparticle)
            do jparticle = istart, istart + nassigned - 1
                cnt = cnt + nassigned_hmc(jparticle)
            end do
        end do
        if (cnt /= nparticle) then
            print *, "wrong: ", theta
            print *, "cnt: ", cnt, "should be the same as ", nparticle
            print *, "sum_score: ", sum_score
            print *, "score_ls: ", score_ls
            print *, "weights_hmc: ", weights_hmc
            print *, "jparticle: "
            do iparticle = 1, nparticle, tuning_factor
                istart = id_start(iparticle)
                nassigned = assigned_num(iparticle)
                do jparticle = istart, istart + nassigned - 1
                    print *, jparticle
                end do
            end do
            print *, "jparticle: end"

            stop
        end if

        cnt = 0
        do iparticle = 1, nparticle, tuning_factor
            istart = id_start(iparticle)
            nassigned = assigned_num(iparticle)
            do jparticle = istart, istart + nassigned - 1
                do kparticle = 1, nassigned_hmc(jparticle)
                    cnt = cnt + 1
                    dtau_ls_new(cnt) = dtau_ls(jparticle)
                    ntau_ls_new(cnt) = ntau_ls(jparticle)
                    ! print *, ntau_ls(jparticle)
                end do
            end do
        end do

!$omp parallel do private(iparticle)
        do iparticle = 1, nparticle
            dtau_ls(iparticle) = dtau_ls_new(iparticle)
            ntau_ls(iparticle) = ntau_ls_new(iparticle)
        end do
!$omp end parallel do

        dtau_mean = 0d0
        ntau_mean = 0
!$omp parallel do private(iparticle) reduction(+:dtau_mean, ntau_mean)
        do iparticle = 1, nparticle
            dtau_mean = dtau_mean + dtau_ls(iparticle)
            ntau_mean = ntau_mean + ntau_ls(iparticle)
        end do
!$omp end parallel do
        dtau_mean = dtau_mean/d_nparticle
        ntau_mean = ntau_mean/nparticle
        log_dtau_upper = log(dtau_mean) + log(5d0)
        log_dtau_lower = log(dtau_mean) - log(5d0)
        ! ntau_upper = min(ntau_mean + 20, 100)
        ! ntau_lower = max(ntau_mean - 20, 1)
        ! print *, "ntau range: ", ntau_lower, ntau_upper
        ! print *, "ntau mean: ", ntau_mean
        ! print *, "dtau range: ", exp(log_dtau_lower), exp(log_dtau_upper)

    end subroutine slip_hmc_tuning

    subroutine slip_mcmc_sampling(gamma, particles, particles_new, &
                                  likelihood_ls, likelihood_ls_new, prior_ls, &
                                  prior_ls_new, cov, assigned_num, id_start, &
                                  nparticle, ndim, dvec, sigma2_full, alpha2_full, &
                                  theta, nplane, nxi, neta, &
                                  gmat, log_sigma_sar2, log_sigma_gnss2, &
                                  ngnss, nobs, ndof, lmat_index, lmat_val, &
                                  nnode, max_slip, st_rand_ls, metropolis_ls, &
                                  particle_cur, particle_cand, st_rand, gsvec, lsvec, &
                                  gsgmat_l, gsdvec, dsd, nsar_total)
        implicit none
        double precision, intent(in) :: max_slip, dvec(:), &
            lmat_val(:, :), gmat(:, :), sigma2_full(:), alpha2_full(:), theta(:), &
            log_sigma_sar2, log_sigma_gnss2, gamma, cov(:, :), gsgmat_l(:, :), &
            gsdvec(:), dsd
        integer, intent(in) :: nnode, ndof, ngnss, nobs, nplane, nxi, neta, &
                               nparticle, ndim, lmat_index(:, :), assigned_num(:), &
                               nsar_total
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
            call random_number_correction(metropolis)
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
                    ! !   non negative constraints
                    ! if (particle_cand(idim) < 0d0) then
                    !     particle_cand(idim) = -particle_cand(idim)
                    ! end if
                    ! !   max slip constraints
                    ! if (particle_cand(idim) > max_slip) then
                    !     particle_cand(idim) = 2*max_slip - particle_cand(idim)
                    ! end if
                end do

                ! calculate negative log likelihood/prior/posterior
                ! of the proposed configuration
                likelihood_cand = slip_calc_likelihood( &
                                  particle_cand, dvec, sigma2_full, gmat, &
                                  log_sigma_sar2, log_sigma_gnss2, ngnss, &
                                  gsvec, nobs, ndof, gsgmat_l, gsdvec, dsd, nsar_total)
                prior_cand = slip_calc_prior(particle_cand, alpha2_full, theta, nplane, &
                                             lmat_index, lmat_val, lsvec, nnode)

                !  metropolis test
                metropolis = metropolis_ls(jparticle)
                if (post_cur - post_cand > log(metropolis)) then
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

    subroutine slip_calc_gsd_gsg(gmat, dvec, gsdvec, gsgmat, &
                                 nobs, ndof, npath, nsar_index, &
                                 nsar_total, sigma_sar_mat, &
                                 sigma_sar2, sigma_gnss2, &
                                 obs_sigma, gmat_arr, dsd)
        implicit none
        integer, intent(in) :: nobs, ndof, npath, nsar_index(:), nsar_total
        type(mat), intent(in) :: sigma_sar_mat(:), gmat_arr(:)
        double precision, intent(in) :: gmat(:, :), dvec(:), &
                                        sigma_sar2, sigma_gnss2, obs_sigma(:)
        double precision, intent(inout) :: gsdvec(:), gsgmat(:, :), dsd
        double precision :: gtsmat(2*ndof, nobs), sdvec(nobs)
        integer :: idim, jdim, iobs, jobs, n, ipath, info
        double precision :: st_time, en_time, dtmp
    
        ! open(10, file="tmp1/sigma1", status="replace")
        ! do iobs = 1, sigma_sar_mat(1)%nrow
        !     do jobs = 1, sigma_sar_mat(1)%nrow
        !         write(10, "(e20.10)", advance="no") sigma_sar_mat(1)%body(iobs, jobs)/sigma_sar2
        !     end do
        !     write(10, *)
        ! end do
        ! close(10)

        ! open(10, file="tmp1/sigma2", status="replace")
        ! do iobs = 1, sigma_sar_mat(2)%nrow
        !     do jobs = 1, sigma_sar_mat(2)%nrow
        !         write(10, "(e20.10)", advance="no") sigma_sar_mat(2)%body(iobs, jobs)/sigma_sar2
        !     end do
        !     write(10, *)
        ! end do
        ! close(10)

        ! open(10, file="tmp1/sigma3", status="replace")
        ! do iobs = 1, sigma_sar_mat(3)%nrow
        !     do jobs = 1, sigma_sar_mat(3)%nrow
        !         write(10, "(e20.10)", advance="no") sigma_sar_mat(3)%body(iobs, jobs)/sigma_sar2
        !     end do
        !     write(10, *)
        ! end do
        ! close(10)

        ! open(10, file="tmp1/sigma4", status="replace")
        ! do iobs = 1, sigma_sar_mat(4)%nrow
        !     do jobs = 1, sigma_sar_mat(4)%nrow
        !         write(10, "(e20.10)", advance="no") sigma_sar_mat(4)%body(iobs, jobs)/sigma_sar2
        !     end do
        !     write(10, *)
        ! end do
        ! close(10)

        ! open(10, file="tmp1/sigmagnss", status="replace")
        ! do iobs = nsar_index(npath + 1),  nobs
        !     do jobs = nsar_index(npath + 1),  nobs
        !         if (iobs == jobs) then
        !             write(10, "(e20.10)", advance="no") 1d0/(obs_sigma(iobs)**2*sigma_gnss2)
        !         else 
        !             write(10, "(e20.10)", advance="no") 0d0
        !         end if
        !     end do
        !     write(10, *)
        ! end do
        ! close(10)

        ! open(10, file="tmp1/gmat", status="replace")
        ! do iobs = 1, nobs
        !     do idim = 1, 2*ndof
        !         write(10, "(e20.10)", advance="no") gmat(iobs, idim)
        !     end do
        !     write(10, *)
        ! end do

        ! open(10, file="tmp1/dvec", status="replace")
        ! do iobs = 1, nobs
        !     write(10, "(e20.10)", advance="no") dvec(iobs)
        ! end do
        ! print *, "sigma_sar2: ", sigma_sar2

        ! sd
        sdvec = 0d0
        ! sar
        do ipath = 1, npath
            n = sigma_sar_mat(ipath)%nrow
            call dgemv("n", n, n, 1d0/sigma_sar2, &
                       sigma_sar_mat(ipath)%body, n, dvec(nsar_index(ipath)), &
                       1, 0d0, sdvec(nsar_index(ipath)), 1)
        end do
        ! gnss
        do iobs = nsar_total + 1, nobs
            sdvec(iobs) = dvec(iobs)/(sigma_gnss2 * obs_sigma(iobs)**2)
        end do

        ! gsd
        gsdvec = 0d0
        do idim = 1, 2*ndof
            do iobs = 1, nobs
                gsdvec(idim) = gsdvec(idim) + gmat(iobs, idim)*sdvec(iobs)
            end do
        end do

        ! dsd
        dsd = 0d0
        do iobs = 1, nobs
            dsd = dsd + dvec(iobs) * sdvec(iobs)
        end do

        gsgmat = 0d0
        gtsmat = 0d0
        ! gts
        ! sar
        do ipath = 1, npath
            call dgemm("t", "n", gmat_arr(ipath)%ncol, gmat_arr(ipath)%nrow, &
                       gmat_arr(ipath)%nrow, 1d0/sigma_sar2, gmat_arr(ipath)%body, &
                       gmat_arr(ipath)%nrow, sigma_sar_mat(ipath)%body, &
                       sigma_sar_mat(ipath)%nrow, 0d0, gtsmat(1, nsar_index(ipath)), &
                       2*ndof)
        end do
        ! gnss
        do iobs = nsar_total + 1, nobs
            do idim = 1, 2*ndof
                gtsmat(idim, iobs) = gmat(iobs, idim)/(obs_sigma(iobs)**2 * sigma_gnss2)
            end do
        end do
        ! gsg
        call dgemm('n', 'n', 2*ndof, 2*ndof, nobs, 1d0, gtsmat, &
                   2*ndof, gmat, nobs, 0d0, gsgmat, 2*ndof)

        ! ! cholesky
        ! gsgmat_l = gsgmat
        ! call dpotrf("u", 2*ndof, gsgmat_l, 2*ndof, info)

        ! open(10, file="tmp1/gsgmat", status="replace")
        ! do idim = 1, 2*ndof
        !     do jdim = 1, 2*ndof
        !         write(10, "(e20.10)", advance="no") gsgmat(idim, jdim)
        !     end do
        !     write(10, *)
        ! end do
        ! close(10)

        ! open(10, file="tmp1/gsgmat_l", status="replace")
        ! do idim = 1, 2*ndof
        !     do jdim = 1, 2*ndof
        !         write(10, "(e20.10)", advance="no") gsgmat_l(idim, jdim)
        !     end do
        !     write(10, *)
        ! end do
        ! close(10)

        ! open(10, file="tmp1/gsdvec", status="replace")
        ! do idim = 1, 2*ndof
        !     write(10, "(e20.10)", advance="no") gsdvec(idim)
        ! end do
        ! close(10)

        ! open(10, file="tmp1/dsd", status="replace")
        !     write(10, "(e20.10)", advance="no") dsd
        ! close(10)
        ! stop
    end subroutine slip_calc_gsd_gsg

    subroutine slip_smc_exec(particles, particles_new, &
                             nparticle, ndim, flag_output, &
                             output_path, llmat, &
                             max_slip, dvec, &
                             sigma2_full, alpha2_full, theta, nplane, &
                             gmat, log_sigma_sar2, log_sigma_gnss2, &
                             ngnss, nobs, ndof, lmat_index, &
                             lmat_val, ltmat_index, ltmat_val, &
                             nnode, likelihood_ls, &
                             prior_ls, weights, mean, &
                             cov, likelihood_ls_new, &
                             prior_ls_new, assigned_num, id_start, &
                             id_dof, st_rand_ls, metropolis_ls, &
                             gsvec, lsvec, particle_cur, particle_cand, &
                             st_rand, neglog_ret, npath, nsar_index, &
                             nsar_total, sigma_sar_mat, sigma_sar2, sigma_gnss2, &
                             obs_sigma, gmat_arr)
        implicit none
        double precision, intent(in) ::  max_slip, dvec(:), &
            lmat_val(:, :), ltmat_val(:, :), llmat(:, :), gmat(:, :), &
            sigma2_full(:), alpha2_full(:), &
            log_sigma_sar2, log_sigma_gnss2, theta(:), &
            sigma_sar2, sigma_gnss2, obs_sigma(:)
        integer, intent(in) :: nplane, nnode, ndof, ngnss, nobs, &
                               nparticle, ndim, flag_output, id_dof(:), &
                               lmat_index(:, :), ltmat_index(:, :), &
                               npath, nsar_total, nsar_index(:)
        character(*), intent(in) :: output_path
        type(mat), intent(in) :: sigma_sar_mat(:), gmat_arr(:)
        double precision, intent(inout) ::  particles(:, :), &
            particles_new(:, :), likelihood_ls(:), &
            prior_ls(:), weights(:), mean(:), cov(:, :), &
            likelihood_ls_new(:), prior_ls_new(:), &
            st_rand_ls(:, :), metropolis_ls(:), gsvec(:), lsvec(:), &
            particle_cur(:), particle_cand(:), st_rand(:), neglog_ret
        integer, intent(inout) ::  assigned_num(:), id_start(:)
        integer :: iparticle, idim, idof, inode, idirection
        integer :: iter, iiter
        double precision :: gamma, neglog_evidence
        double precision, allocatable :: slip(:, :)
        double precision :: st_time1, st_time2, en_time1, en_time2
        character(len=200) :: iter_char, filename
        character(len=200) :: output_dir
        ! HMC
        double precision :: cov_diag(ndim)
        ! tuning HMC
        double precision :: log_dtau_upper, log_dtau_lower, dtau_ls(nparticle)
        integer :: ntau_upper, ntau_lower, ntau_ls(nparticle)

        double precision :: gsdvec(ndof*2), gsgmat(ndof*2, ndof*2), dsd
        integer :: tuning_factor
        double precision :: metropolis, v1, v2

        output_dir = "./output_slip/"
        ! log_dtau_upper = 1d1
        log_dtau_upper = 0d0
        log_dtau_lower = -1d1
        ntau_lower = 1
        ntau_upper = 100
        tuning_factor = 5
        iter = 0

        call slip_calc_gsd_gsg(gmat, dvec, gsdvec, gsgmat, nobs, ndof, npath, &
                               nsar_index, nsar_total, sigma_sar_mat, &
                               sigma_sar2, sigma_gnss2, obs_sigma, gmat_arr, dsd)

        ! sampling from the prior distribution
        call slip_gen_init_particles(particles, likelihood_ls, prior_ls, nparticle, &
                                     ndim, llmat, max_slip, dvec, &
                                     sigma2_full, alpha2_full, theta, nplane, &
                                     gmat, log_sigma_sar2, &
                                     log_sigma_gnss2, ngnss, nobs, ndof, &
                                     lmat_index, lmat_val, nnode, particle_cur, gsvec, lsvec, &
                                     gsdvec, dsd, nsar_total, gsgmat)
        ! ! output result of stage 0(disabled)
        ! write (iter_char, "(i0)") iter
        ! filename = trim(trim(output_dir)//trim(iter_char)//".csv")
        ! open (17, file=filename, status='replace')
        ! do iparticle = 1, nparticle
        !     do idim = 1, ndim
        !         write (17, "(f12.5)", advance="no") particles(idim, iparticle)
        !     end do
        !     write (17, "(f12.5)", advance="no") likelihood_ls(iparticle)
        !     write (17, "(f12.5)", advance="no") prior_ls(iparticle)
        !     write (17, *)
        ! end do
        ! close (17)
        iter = iter + 1

        gamma = 0d0
        ! product of S_j (sum for negative log value)
        neglog_ret = 0d0

        do while (1d0 - gamma > 10d-8)
            ! print *, "iter: ", iter
            ! do iiter = 1, 10
            st_time1 = omp_get_wtime()
            ! S_j

            ! find the gamma such that c.o.v of weights = 0.5
            gamma = slip_find_next_gamma(gamma, likelihood_ls, weights, &
                                         neglog_evidence, nparticle)
            print *, "gamma: ", gamma
            print *, "evidence: ", neglog_evidence
            neglog_ret = neglog_ret + neglog_evidence
            if (iter > 200) then
                neglog_ret = 1d10
                exit
            end if

            ! normalize weights(sum of weights needs to be 1)
            call slip_normalize_weights(weights, nparticle)

            ! calculate mean and covariance of the samples
            call slip_calc_mean_particles(particles, weights, &
                                          mean, nparticle, ndim)
            call slip_calc_cov_particles(particles, weights, mean, cov, &
                                         nparticle, ndim, cov_diag)
            call slip_resample_particles(nparticle, weights, assigned_num)

            st_time2 = omp_get_wtime()

            call slip_hmc_tuning(gamma, particles, particles_new, &
                                 likelihood_ls, likelihood_ls_new, prior_ls, &
                                 prior_ls_new, cov, cov_diag, assigned_num, id_start, &
                                 nparticle, ndim, dvec, sigma2_full, alpha2_full, &
                                 theta, nplane, &
                                 gmat, log_sigma_sar2, log_sigma_gnss2, &
                                 ngnss, nobs, ndof, lmat_index, lmat_val, &
                                 ltmat_index, ltmat_val, &
                                 nnode, max_slip, st_rand_ls, metropolis_ls, &
                                 particle_cur, particle_cand, st_rand, gsvec, lsvec, &
                                 log_dtau_upper, log_dtau_lower, ntau_upper, ntau_lower, &
                                 dtau_ls, ntau_ls, &
                                 gsdvec, gsgmat, tuning_factor, dsd, nsar_total)

            call slip_hmc_sampling(gamma, particles, particles_new, &
                                   likelihood_ls, likelihood_ls_new, prior_ls, &
                                   prior_ls_new, cov, cov_diag, assigned_num, id_start, &
                                   nparticle, ndim, dvec, sigma2_full, alpha2_full, &
                                   theta, nplane, &
                                   gmat, log_sigma_sar2, log_sigma_gnss2, &
                                   ngnss, nobs, ndof, lmat_index, lmat_val, &
                                   ltmat_index, ltmat_val, &
                                   nnode, max_slip, st_rand_ls, metropolis_ls, &
                                   particle_cur, particle_cand, st_rand, gsvec, lsvec, &
                                   dtau_ls, ntau_ls, gsdvec, gsgmat, ntau_upper, dsd, nsar_total)

            ! call slip_hmc_sampling(gamma, particles, particles_new, &
            !                        likelihood_ls, likelihood_ls_new, prior_ls, &
            !                        prior_ls_new, cov, cov_diag, assigned_num, id_start, &
            !                        nparticle, ndim, dvec, sigma2_full, alpha2_full, &
            !                        theta, nplane, nxi, neta, &
            !                        gmat, log_sigma_sar2, log_sigma_gnss2, nsar, &
            !                        ngnss, nobs, ndof, lmat_index, lmat_val, &
            !                        ltmat_index, ltmat_val, &
            !                        nnode, max_slip, st_rand_ls, metropolis_ls, &
            !                        particle_cur, particle_cand, st_rand, gsvec, lsvec, &
            !                        dtau_ls, ntau_ls, gsdvec, gsgmat, ntau_upper)

            ! call slip_mcmc_sampling(gamma, particles, particles_new, &
            !                         likelihood_ls, likelihood_ls_new, &
            !                         prior_ls, prior_ls_new, &
            !                         cov, assigned_num, id_start, nparticle, &
            !                         ndim, dvec, sigma2_full, alpha2_full, &
            !                         theta, nplane, nxi, neta, gmat, log_sigma_sar2, &
            !                         log_sigma_gnss2, nsar, ngnss, nobs, ndof, &
            !                         lmat_index, lmat_val, nnode, max_slip, &
            !                         st_rand_ls, metropolis_ls, particle_cur, &
            !                         particle_cand, st_rand, gsvec, lsvec)
            en_time1 = omp_get_wtime()
            ! print *, "loop total: ", en_time1 - st_time1
            ! print *, "hmc: ", en_time1 - st_time2
            ! ! output result of stage iter (disabled)
            ! write (iter_char, "(i0)") iter
            ! filename = trim(trim(output_dir)//trim(iter_char)//".csv")
            ! open (17, file=filename, status='replace')
            ! do iparticle = 1, nparticle
            !     do idim = 1, ndim
            !         write (17, "(f12.5)", advance="no") particles(idim, iparticle)
            !     end do
            !     write (17, "(f12.5)", advance="no") likelihood_ls(iparticle)
            !     write (17, "(f12.5)", advance="no") prior_ls(iparticle)
            !     write (17, *)
            ! end do
            ! close (17)
            iter = iter + 1
        end do
        ! print *, "max iter:", iter
        ! print *, "neglog", neglog_ret

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
                ! write (17, "(f12.5)") likelihood_ls(iparticle) + prior_ls(iparticle)
                write (17, "(f12.5)", advance="no") likelihood_ls(iparticle) 
                write (17, "(f12.5)") prior_ls(iparticle)
            end do
            close (17)
            deallocate (slip)
        end if
        ! print *, "iter: ", iter, " neglog_ret: ", neglog_ret

    end subroutine
end module smc_slip
