module smc_fault
    use init
    use gfunc
    use smc_slip
    use mpi
    use omp_lib
    implicit none
contains
    subroutine fault_sample_init_particles(particles, nparticle, ndim, range)
        implicit none

        double precision, intent(inout) :: particles(:, :)
        integer, intent(in) :: nparticle, ndim
        double precision, intent(in) :: range(:, :)
        integer :: idim, iparticle
        double precision ::  x, xmin, xmax, randr
        do iparticle = 1, nparticle
            do idim = 1, ndim
                xmin = range(1, idim)
                xmax = range(2, idim)
                call random_number(randr)
                x = xmin + randr*(xmax - xmin)
                particles(idim, iparticle) = x
            end do
        end do
    end subroutine fault_sample_init_particles

    double precision function fault_calc_likelihood( &
        particle, nxi, neta, nnode, ndof, nsar, &
        ngnss, nobs, cny_fault, coor_fault, &
        node_to_elem_val, node_to_elem_size, id_dof, luni, &
        lmat, lmat_index, lmat_val, llmat, &
        gmat, slip_dist, obs_points, &
        obs_unitvec, obs_sigma, sigma2_full, &
        target_id_val, node_id_in_patch, xinode, etanode, &
        uxinode, uetanode, r1vec, r2vec, &
        nvec, response_dist, uobs, uret, &
        slip_particles, slip_particles_new, &
        nparticle_slip, max_slip, dvec, &
        slip_likelihood_ls, slip_prior_ls, slip_weights, &
        slip_mean, slip_cov, slip_likelihood_ls_new, &
        slip_prior_ls_new, slip_assigned_num, slip_id_start, &
        slip_st_rand_ls, slip_metropolis_ls, gsvec, lsvec, &
        slip_particle_cur, slip_particle_cand, slip_st_rand, &
        flag_output, output_path)
        implicit none
        double precision, intent(in) :: particle(:), obs_points(:, :), &
            obs_unitvec(:, :), obs_sigma(:), max_slip, dvec(:)
        integer, intent(in) :: nxi, neta, nnode, ndof, nsar, ngnss, nobs, &
                               nparticle_slip, flag_output
        character(*), intent(in) :: output_path
        double precision, intent(inout) :: coor_fault(:, :), luni(:, :), &
            lmat(:, :), lmat_val(:, :), llmat(:, :), gmat(:, :), slip_dist(:, :), &
            sigma2_full(:), xinode(:), etanode(:), uxinode(:), uetanode(:), &
            r1vec(:), r2vec(:), nvec(:), response_dist(:, :), uobs(:), uret(:), &
            slip_particles(:, :), slip_particles_new(:, :), slip_likelihood_ls(:), &
            slip_prior_ls(:), slip_weights(:), slip_mean(:), slip_cov(:, :), &
            slip_likelihood_ls_new(:), slip_prior_ls_new(:), &
            slip_st_rand_ls(:, :), slip_metropolis_ls(:), gsvec(:), lsvec(:), &
            slip_particle_cur(:), slip_particle_cand(:), slip_st_rand(:)
        integer, intent(inout) :: cny_fault(:, :), node_to_elem_val(:, :), &
                                  node_to_elem_size(:), id_dof(:), &
                                  lmat_index(:, :), target_id_val(:), &
                                  node_id_in_patch(:), slip_assigned_num(:), &
                                  slip_id_start(:)
        integer :: i, j
        double precision :: xf, yf, zf, strike, dip, log_sigma_sar2, &
            log_sigma_gnss2, log_alpha2, lxi, leta, dxi, deta, neglog_sum
!     // double st_time, en_time
!     // st_time = MPI_Wtime()
        xf = particle(1)
        yf = particle(2)
        zf = particle(3)
        strike = particle(4)
        dip = particle(5)
        log_sigma_sar2 = particle(6)
        log_sigma_gnss2 = particle(7)
        log_alpha2 = particle(8)
        lxi = particle(9)
        leta = particle(10)

        ! set fault geometry
        call discretize_fault(lxi, leta, nxi, neta, cny_fault, coor_fault, &
                              node_to_elem_val, node_to_elem_size, id_dof)
        dxi = lxi/nxi
        deta = leta/neta

        ! calculate laplacian matrix L
        call gen_laplacian(nnode, nxi, neta, dxi, deta, id_dof, ndof, luni, lmat)
        ! sparse matrix form of lmat
        call gen_sparse_lmat(lmat, lmat_index, lmat_val, nnode, ndof)
        ! matrix L^T L
        call calc_ll(llmat, lmat, nnode, ndof)

        ! Calculate greens function for the sampled fault
        call calc_greens_func(gmat, slip_dist, cny_fault, coor_fault, obs_points, &
                              obs_unitvec, leta, xf, yf, zf, strike, dip, &
                              node_to_elem_val, node_to_elem_size, id_dof, nsar, ngnss, &
                              nobs, nnode, ndof, target_id_val, node_id_in_patch, xinode, &
                              etanode, uxinode, uetanode, r1vec, r2vec, nvec, &
                              response_dist, uobs, uret)

        ! diag component of sigma
        ! (variance matrix for the likelihood function of slip)
        do i = 1, nsar
            sigma2_full(i) = obs_sigma(i)**2d0*exp(log_sigma_sar2)
        end do
        do i = 1, ngnss
            do j = 1, 3
                sigma2_full(nsar + 3*(i - 1) + j) = &
                    obs_sigma(nsar + 3*(i - 1) + j)**2d0*exp(log_sigma_gnss2)
            end do
        end do

!     // en_time = MPI_Wtime()
!     // printf("other: %f\n", en_time - st_time)

!     // st_time = MPI_Wtime()
        ! Sequential Monte Carlo sampling for slip
        ! calculate negative log of likelihood
        call slip_smc_exec( &
            slip_particles, slip_particles_new, nparticle_slip, 2*ndof, &
            flag_output, output_path, llmat, log_alpha2, max_slip, dvec, sigma2_full, &
            gmat, log_sigma_sar2, log_sigma_gnss2, nsar, ngnss, nobs, ndof, &
            lmat_index, lmat_val, nnode, slip_likelihood_ls, slip_prior_ls, &
            slip_weights, slip_mean, slip_cov, slip_likelihood_ls_new, &
            slip_prior_ls_new, slip_assigned_num, slip_id_start, id_dof, &
            slip_st_rand_ls, slip_metropolis_ls, gsvec, lsvec, slip_particle_cur, &
            slip_particle_cand, slip_st_rand, neglog_sum)
        fault_calc_likelihood = neglog_sum
    end function fault_calc_likelihood

    subroutine work_eval_init_particles(work_size, ndim, particle_cur, &
                                        work_particles, work_likelihood_ls, nxi, neta, nnode, ndof, nsar, ngnss, nobs, &
                                        cny_fault, coor_fault, node_to_elem_val, node_to_elem_size, &
                                        id_dof, luni, lmat, lmat_index, lmat_val, llmat, gmat, &
                                        slip_dist, obs_points, obs_unitvec, obs_sigma, sigma2_full, &
                                        target_id_val, node_id_in_patch, xinode, etanode, uxinode, &
                                        uetanode, r1vec, r2vec, nvec, response_dist, uobs, uret, &
                                        slip_particles, slip_particles_new, nparticle_slip, max_slip, &
                                        dvec, slip_likelihood_ls, slip_prior_ls, slip_weights, slip_mean, &
                                        slip_cov, slip_likelihood_ls_new, slip_prior_ls_new, slip_assigned_num, &
                                        slip_id_start, slip_st_rand_ls, slip_metropolis_ls, gsvec, lsvec, &
                                        slip_particle_cur, slip_particle_cand, slip_st_rand)
        implicit none
        integer, intent(in) :: work_size, ndim, nxi, neta, nnode, &
                               ndof, nsar, ngnss, nobs, nparticle_slip
        double precision, intent(in) :: work_particles(:, :), obs_points(:, :), &
            obs_unitvec(:, :), obs_sigma(:), max_slip, dvec(:)
        integer, intent(inout) :: cny_fault(:, :), node_to_elem_val(:, :), node_to_elem_size(:), &
                                  id_dof(:), lmat_index(:, :), target_id_val(:), node_id_in_patch(:), &
                                  slip_assigned_num(:), slip_id_start(:)
        double precision, intent(inout) :: work_likelihood_ls(:), particle_cur(:), coor_fault(:, :), luni(:, :), lmat(:, :), &
            lmat_val(:, :), llmat(:, :), gmat(:, :), slip_dist(:, :), sigma2_full(:), xinode(:), &
            etanode(:), uxinode(:), uetanode(:), r1vec(:), r2vec(:), nvec(:), response_dist(:, :), &
            uobs(:), uret(:), slip_particles(:, :), slip_particles_new(:, :), &
            slip_likelihood_ls(:), slip_prior_ls(:), slip_weights(:), slip_mean(:), slip_cov(:, :), &
            slip_likelihood_ls_new(:), slip_prior_ls_new(:), slip_st_rand_ls(:, :), slip_metropolis_ls(:), &
            gsvec(:), lsvec(:), slip_particle_cur(:), slip_particle_cand(:), slip_st_rand(:)
        integer :: iparticle, idim
        double precision :: likelihood
        do iparticle = 1, work_size
            do idim = 1, ndim
                particle_cur(idim) = work_particles(idim, iparticle)
            end do
            ! calculate negative log likelihood for the sample
            likelihood = fault_calc_likelihood( &
                         particle_cur, nxi, neta, nnode, ndof, nsar, ngnss, nobs, cny_fault, &
                         coor_fault, node_to_elem_val, node_to_elem_size, id_dof, luni, lmat, &
                         lmat_index, lmat_val, llmat, gmat, slip_dist, obs_points, &
                         obs_unitvec, obs_sigma, sigma2_full, target_id_val, node_id_in_patch, &
                         xinode, etanode, uxinode, uetanode, r1vec, r2vec, nvec, response_dist, &
                         uobs, uret, slip_particles, slip_particles_new, &
                         nparticle_slip, max_slip, dvec, slip_likelihood_ls, slip_prior_ls, &
                         slip_weights, slip_mean, slip_cov, slip_likelihood_ls_new, &
                         slip_prior_ls_new, slip_assigned_num, slip_id_start, slip_st_rand_ls, &
                         slip_metropolis_ls, gsvec, lsvec, slip_particle_cur, &
                         slip_particle_cand, slip_st_rand, 0, "")
            work_likelihood_ls(iparticle) = likelihood
        end do
    end subroutine work_eval_init_particles

    subroutine fault_box_muller(p, q)
        implicit none

        double precision, intent(inout) :: p, q
        double precision :: pi, r, s
        pi = 2d0*asin(1d0)
        ! uniform random numbers between 0 and 1
        call random_number(r)
        call random_number(s)
        ! Gaussian random numbers,
        ! with weights proportional to eˆ{-pˆ2/2}and eˆ{-qˆ2/2}
        p = sqrt(-2d0*log(r))*sin(2d0*pi*s)
        q = sqrt(-2d0*log(r))*cos(2d0*pi*s)
    end subroutine fault_box_muller

    subroutine calc_mean_std_vector(vec, size, mean, std)
        implicit none

        double precision, intent(in) :: vec(:)
        integer, intent(in) :: size
        double precision, intent(inout) :: mean, std
        integer i
        mean = 0d0
        do i = 1, size
            mean = mean + vec(i)
        end do
        mean = mean/size
        std = 0d0
        do i = 1, size
            std = std + (vec(i) - mean)**2d0
        end do
        std = std/size
        std = sqrt(std)
    end subroutine calc_mean_std_vector

    subroutine fault_find_next_gamma(gamma, gamma_prev, likelihood_ls, &
                                     weights, nparticle)
        implicit none

        double precision, intent(inout) :: gamma
        double precision, intent(in) :: gamma_prev
        double precision, intent(in) :: likelihood_ls(:)
        double precision, intent(inout) :: weights(:)
        integer, intent(in) :: nparticle
        integer :: iparticle
        double precision :: min_likelihood, likelihood, cv_threshold
        double precision :: lower, upper, err
        double precision :: diff_gamma, mean, std, cv
        cv_threshold = 1d0
        ! find minimum of negative log likelihood
        min_likelihood = likelihood_ls(1)
        do iparticle = 1, nparticle
            min_likelihood = min(min_likelihood, likelihood_ls(iparticle))
        end do
        ! binary search for the next gamma
        ! such that c.o.v of the weight is equivalent to cv_threashold
        lower = gamma_prev
        upper = 1d0
        err = 1d0
        do while (err > 10d0**(-8d0))
            gamma = (lower + upper)/2d0
            diff_gamma = gamma - gamma_prev
            do iparticle = 1, nparticle
                likelihood = likelihood_ls(iparticle)
                ! extract min_likelihood to avoid underflow of the weight
                weights(iparticle) = exp(-diff_gamma*(likelihood - min_likelihood))
            end do

            ! calculate c.o.v = mean/std
            call calc_mean_std_vector(weights, nparticle, mean, std)
            cv = std/mean
            if (cv > cv_threshold) then
                upper = gamma
            else
                lower = gamma
            end if
            err = abs(cv - cv_threshold)
            if (abs(gamma - 1) < 10d0**(-8d0)) then
                exit
            end if
        end do
    end subroutine fault_find_next_gamma

    subroutine fault_normalize_weights(weights, nparticle)
        implicit none

        double precision, intent(inout) :: weights(:)
        integer, intent(in) :: nparticle
        double precision :: sum
        integer :: iparticle
        sum = 0d0
        do iparticle = 1, nparticle
            sum = sum + weights(iparticle)
        end do
        do iparticle = 1, nparticle
            weights(iparticle) = weights(iparticle)/sum
        end do
    end subroutine fault_normalize_weights

    subroutine fault_calc_mean_particles(particles, weights, mean, &
                                         nparticle, ndim)
        implicit none

        double precision, intent(in) :: particles(:, :), weights(:)
        double precision, intent(inout) :: mean(:)
        integer, intent(in) :: nparticle, ndim
        integer :: iparticle, idim
        do idim = 1, ndim
            mean(idim) = 0d0
        end do
        do iparticle = 1, nparticle
            do idim = 1, ndim
                mean(idim) = mean(idim) + &
                             weights(iparticle)*particles(idim, iparticle)
            end do
        end do
    end subroutine fault_calc_mean_particles

    subroutine fault_calc_cov_particles(particles, weights, mean, &
                                        cov, nparticle, ndim)
        implicit none

        double precision, intent(in) :: particles(:, :), weights(:), mean(:)
        double precision, intent(inout) :: cov(:, :)
        integer, intent(in) :: nparticle, ndim
        integer :: iparticle, idim, jdim, ierr
        double precision :: weight

        do jdim = 1, ndim
            do idim = 1, ndim
                cov(idim, jdim) = 0d0
            end do
        end do

        do iparticle = 1, nparticle
            weight = weights(iparticle)
            do jdim = 1, ndim
                do idim = jdim, ndim
                    cov(idim, jdim) = cov(idim, jdim) + &
                                      weight* &
                                      (particles(idim, iparticle) - mean(idim))* &
                                      (particles(jdim, iparticle) - mean(jdim))
                end do
            end do
        end do

        write (*, "(a)", advance="no") "mean: "
        do idim = 1, ndim
            write (*, "(f12.5)", advance="no") mean(idim)
        end do
        write (*, *)

        write (*, "(a)", advance="no") "cov: "
        do idim = 1, ndim
            write (*, "(f12.5)", advance="no") cov(idim, idim)
        end do
        write (*, *)

        do jdim = 1, ndim
            do idim = 1, ndim
                cov(idim, jdim) = cov(idim, jdim)*0.04d0
            end do
        end do

        ! LAPACK function for LU decomposition of matrix
        call dpotrf('L', ndim, cov, ndim, ierr)
    end subroutine fault_calc_cov_particles

    subroutine fault_resample_particles(nparticle, weights, assigned_num)
        implicit none

        integer, intent(in) :: nparticle
        double precision, intent(in) :: weights(:)
        integer, intent(inout) :: assigned_num(:)
        integer :: iparticle
        double precision ::  d_nparticle, u
        d_nparticle = nparticle
        ! systematic residual resampling
        call random_number(u)
        u = u/d_nparticle
        ! for(iparticle=0 iparticle < nparticle iparticle + +) {
        do iparticle = 1, nparticle
            assigned_num(iparticle) = &
                floor((weights(iparticle) - u)*d_nparticle) + 1
            u = u + assigned_num(iparticle)/d_nparticle - weights(iparticle)
        end do
    end subroutine fault_resample_particles

    subroutine fault_reorder_to_send(assigned_num, particles, &
                                     tmp_assigned_num, tmp_particles, &
                                     sorted_idx, nparticle, ndim, numprocs, work_size)
        implicit none
        integer, intent(inout) :: assigned_num(:), tmp_assigned_num(:)
        double precision, intent(inout) :: particles(:, :), tmp_particles(:, :)
        integer, intent(inout) :: sorted_idx(:)
        integer, intent(in) :: nparticle, ndim, numprocs, work_size
        integer :: iparticle, idim, iproc
        integer :: id_org, id_new
        do iparticle = 1, nparticle
            tmp_assigned_num(iparticle) = assigned_num(iparticle)
            do idim = 1, ndim
                tmp_particles(idim, iparticle) = &
                    particles(idim, iparticle)
            end do
        end do

        do iparticle = 1, nparticle
            sorted_idx(iparticle) = iparticle
        end do
        call qsort(assigned_num, sorted_idx, nparticle)
        do iparticle = 1, nparticle
            print *, assigned_num(sorted_idx(iparticle))
        end do
        do iproc = 1, numprocs
            do iparticle = 1, work_size
                id_org = sorted_idx((iparticle - 1)*numprocs + iproc)
                id_new = (iproc - 1)*work_size + iparticle
                assigned_num(id_new) = tmp_assigned_num(id_org)
                do idim = 1, ndim
                    particles(idim, id_new) = &
                        tmp_particles(idim, id_org)
                end do
            end do
        end do
    end subroutine fault_reorder_to_send

    subroutine work_mcmc_sampling(work_assigned_num, work_particles, &
                                  work_particles_new, &
                                  work_likelihood_ls, id_start, &
                                  particle_cur, particle_cand, &
                                  st_rand, work_size, ndim, &
                                  cov, gamma, myid, nxi, neta, nnode, ndof, &
                                  nsar, ngnss, nobs, cny_fault, coor_fault, &
                                  node_to_elem_val, node_to_elem_size, id_dof, &
                                  luni, lmat, lmat_index, lmat_val, llmat, gmat, &
                                  slip_dist, obs_points, obs_unitvec, obs_sigma, &
                                  sigma2_full, target_id_val, node_id_in_patch, &
                                  xinode, etanode, uxinode, uetanode, r1vec, r2vec, &
                                  nvec, response_dist, uobs, uret, slip_particles, &
                                  slip_particles_new, nparticle_slip, max_slip, dvec, &
                                  slip_likelihood_ls, slip_prior_ls, slip_weights, &
                                  slip_mean, slip_cov, slip_likelihood_ls_new, slip_prior_ls_new, &
                                  slip_assigned_num, slip_id_start, slip_st_rand_ls, &
                                  slip_metropolis_ls, gsvec, lsvec, slip_particle_cur, &
                                  slip_particle_cand, slip_st_rand)
        implicit none
        integer, intent(in) :: work_assigned_num(:), id_start(:), work_size, ndim, myid, &
                               nxi, neta, nnode, ndof, nsar, ngnss, nobs, nparticle_slip
        double precision, intent(in) :: work_particles(:, :), cov(:, :), gamma, obs_points(:, :), &
            obs_unitvec(:, :), obs_sigma(:), max_slip, dvec(:)
        integer, intent(inout) :: cny_fault(:, :), node_to_elem_val(:, :), node_to_elem_size(:), &
                                  id_dof(:), lmat_index(:, :), target_id_val(:), node_id_in_patch(:), &
                                  slip_assigned_num(:), slip_id_start(:)
        double precision, intent(inout) :: work_particles_new(:, :), work_likelihood_ls(:), &
            particle_cur(:), particle_cand(:), st_rand(:), coor_fault(:, :), luni(:, :), &
            lmat(:, :), lmat_val(:, :), llmat(:, :), gmat(:, :), slip_dist(:, :), sigma2_full(:), &
            xinode(:), etanode(:), uxinode(:), uetanode(:), r1vec(:), r2vec(:), nvec(:), response_dist(:, :), &
            uobs(:), uret(:), slip_particles(:, :), slip_particles_new(:, :), slip_likelihood_ls(:), &
            slip_prior_ls(:), slip_weights(:), slip_mean(:), slip_cov(:, :), slip_likelihood_ls_new(:), &
            slip_prior_ls_new(:), slip_st_rand_ls(:, :), slip_metropolis_ls(:), gsvec(:), lsvec(:), &
            slip_particle_cur(:), slip_particle_cand(:), slip_st_rand(:)
        integer :: iparticle, jparticle, idim, jdim, nassigned, istart
        double precision :: likelihood_cur, likelihood_cand, metropolis
        double precision :: v1, v2

        do iparticle = 1, work_size
            nassigned = work_assigned_num(iparticle)
            istart = id_start(iparticle)
            do idim = 1, ndim
                particle_cur(idim) = work_particles(idim, iparticle)
            end do
            do jparticle = istart, istart + nassigned - 1
                likelihood_cur = fault_calc_likelihood( &
                                 particle_cur, nxi, neta, nnode, ndof, nsar, ngnss, nobs, cny_fault, &
                                 coor_fault, node_to_elem_val, node_to_elem_size, id_dof, luni, lmat, &
                                 lmat_index, lmat_val, llmat, gmat, slip_dist, obs_points, &
                                 obs_unitvec, obs_sigma, sigma2_full, target_id_val, node_id_in_patch, &
                                 xinode, etanode, uxinode, uetanode, r1vec, r2vec, nvec, response_dist, &
                                 uobs, uret, slip_particles, slip_particles_new, &
                                 nparticle_slip, max_slip, dvec, slip_likelihood_ls, slip_prior_ls, &
                                 slip_weights, slip_mean, slip_cov, slip_likelihood_ls_new, &
                                 slip_prior_ls_new, slip_assigned_num, slip_id_start, slip_st_rand_ls, &
                                 slip_metropolis_ls, gsvec, lsvec, slip_particle_cur, &
                                 slip_particle_cand, slip_st_rand, 0, "")
                ! propose particle_cand
                do idim = 1, ndim
                    call fault_box_muller(v1, v2)
                    st_rand(idim) = v1
                end do
                call dtrmv('l', 'n', 'n', &
                           ndim, cov, ndim, st_rand, 1)
                do idim = 1, ndim
                    particle_cand(idim) = particle_cur(idim) + st_rand(idim)
                end do
                ! calculate negative log likelihood of the proposed configuration
                likelihood_cand = fault_calc_likelihood( &
                                  particle_cand, nxi, neta, nnode, ndof, nsar, ngnss, nobs, cny_fault, &
                                  coor_fault, node_to_elem_val, node_to_elem_size, id_dof, luni, lmat, &
                                  lmat_index, lmat_val, llmat, gmat, slip_dist, obs_points, &
                                  obs_unitvec, obs_sigma, sigma2_full, target_id_val, node_id_in_patch, &
                                  xinode, etanode, uxinode, uetanode, r1vec, r2vec, nvec, response_dist, &
                                  uobs, uret, slip_particles, slip_particles_new, &
                                  nparticle_slip, max_slip, dvec, slip_likelihood_ls, slip_prior_ls, &
                                  slip_weights, slip_mean, slip_cov, slip_likelihood_ls_new, &
                                  slip_prior_ls_new, slip_assigned_num, slip_id_start, slip_st_rand_ls, &
                                  slip_metropolis_ls, gsvec, lsvec, slip_particle_cur, &
                                  slip_particle_cand, slip_st_rand, 0, "")
                ! metropolis test and check domain of definition
                call random_number(metropolis)
                if (exp(gamma*(likelihood_cur - likelihood_cand)) > metropolis) then
                    do idim = 1, ndim
                        particle_cur(idim) = particle_cand(idim)
                    end do
                    likelihood_cur = likelihood_cand
                else
                end if

                ! save to new particle list
                work_likelihood_ls(jparticle) = likelihood_cur
                do idim = 1, ndim
                    work_particles_new(idim, jparticle) = &
                        particle_cur(idim)
                end do
            end do
        end do

    end subroutine work_mcmc_sampling

    subroutine fault_smc_exec( &
        output_dir, range, nparticle, ndim, &
        myid, numprocs, nxi, neta, nnode, ndof, nsar, ngnss, nobs, cny_fault, &
        coor_fault, node_to_elem_val, node_to_elem_size, id_dof, luni, lmat, &
        lmat_index, lmat_val, llmat, gmat, slip_dist, obs_points, &
        obs_unitvec, obs_sigma, sigma2_full, target_id_val, node_id_in_patch, &
        xinode, etanode, uxinode, uetanode, r1vec, r2vec, nvec, &
        response_dist, uobs, uret, slip_particles, &
        slip_particles_new, nparticle_slip, max_slip, dvec, gsvec, &
        lsvec, slip_likelihood_ls, slip_prior_ls, slip_weights, slip_mean, &
        slip_cov, slip_likelihood_ls_new, slip_prior_ls_new, &
        slip_st_rand, slip_particle_cur, slip_particle_cand, &
        slip_assigned_num, slip_id_start, slip_st_rand_ls, slip_metropolis_ls)
        implicit none
        double precision, intent(in) :: obs_points(:, :), &
            obs_unitvec(:, :), obs_sigma(:), max_slip, dvec(:), range(:, :)
        integer, intent(in) :: nxi, neta, nnode, ndof, nsar, ngnss, nobs, &
                               nparticle_slip, nparticle, ndim, &
                               myid, numprocs
        double precision, intent(inout) :: coor_fault(:, :), luni(:, :), &
            lmat(:, :), lmat_val(:, :), llmat(:, :), gmat(:, :), slip_dist(:, :), &
            sigma2_full(:), xinode(:), etanode(:), uxinode(:), uetanode(:), &
            r1vec(:), r2vec(:), nvec(:), response_dist(:, :), uobs(:), uret(:), &
            slip_particles(:, :), slip_particles_new(:, :), slip_likelihood_ls(:), &
            slip_prior_ls(:), slip_weights(:), slip_mean(:), slip_cov(:, :), &
            slip_likelihood_ls_new(:), slip_prior_ls_new(:), &
            slip_st_rand_ls(:, :), slip_metropolis_ls(:), gsvec(:), lsvec(:), &
            slip_particle_cur(:), slip_particle_cand(:), slip_st_rand(:)
        integer, intent(inout) :: cny_fault(:, :), node_to_elem_val(:, :), &
                                  node_to_elem_size(:), id_dof(:), &
                                  lmat_index(:, :), target_id_val(:), &
                                  node_id_in_patch(:), slip_assigned_num(:), &
                                  slip_id_start(:)
        character(*), intent(in) :: output_dir
        integer :: work_size, iparticle, iproc, idim, sum_assigned, cnt, iter
        character(len=100) :: iter_char, filename
        double precision :: gamma, gamma_prev
        double precision, allocatable :: particles(:, :)
        ! array for only the master process
        integer, allocatable :: sum_assigned_ls(:)
        double precision, allocatable :: buf_likelihood(:, :), buf_particles_flat(:, :)
        double precision, allocatable :: weights(:), mean(:)
        double precision, allocatable :: likelihood_ls(:)
        integer, allocatable :: assigned_num(:)
        integer, allocatable :: sorted_idx(:)
        double precision, allocatable :: tmp_particles(:, :)
        integer, allocatable :: tmp_assigned_num(:)
        ! for MPI
        double precision, allocatable :: work_particles(:, :), work_particles_new(:, :), &
            work_likelihood_ls(:)
        integer, allocatable :: work_assigned_num(:)

        double precision, allocatable :: particle_cur(:), particle_cand(:)
        integer, allocatable :: id_start(:)
        double precision, allocatable :: cov(:, :), st_rand(:)
        integer :: ierr
        integer, allocatable :: istatus(:)

        ! measuring time
        double precision :: st_time, en_time

        allocate (istatus(mpi_status_size))
        work_size = nparticle/numprocs
        if (myid == 0) then
            allocate (particles(ndim, nparticle))
            allocate (sum_assigned_ls(numprocs))
            allocate (buf_likelihood(50*work_size, numprocs))
            allocate (buf_particles_flat(50*work_size*ndim, numprocs))
            allocate (weights(nparticle))
            allocate (mean(ndim))
            allocate (likelihood_ls(nparticle))
            allocate (assigned_num(nparticle))
            allocate (sorted_idx(nparticle))
            allocate (tmp_particles(nparticle, nparticle))
            allocate (tmp_assigned_num(nparticle))
        else
            allocate (particles(1, 1))
            allocate (likelihood_ls(1))
            allocate (assigned_num(1))
        end if

        allocate (work_particles(ndim, work_size))
        allocate (work_particles_new(ndim, 50*work_size))
        allocate (work_likelihood_ls(50*work_size))
        allocate (work_assigned_num(work_size))
        allocate (particle_cur(ndim))
        allocate (particle_cand(ndim))
        allocate (id_start(work_size))
        allocate (st_rand(ndim))
        allocate (cov(ndim, ndim))
        iter = 0
        if (myid == 0) then
            ! sampling from the prior distribution
            call fault_sample_init_particles(particles, nparticle, ndim, range)
        end if
        call mpi_scatter(particles, ndim*work_size, mpi_double_precision, &
                         work_particles, ndim*work_size, mpi_double_precision, &
                         0, mpi_comm_world, ierr)
        call work_eval_init_particles(work_size, ndim, particle_cur, &
                                      work_particles, work_likelihood_ls, nxi, neta, nnode, ndof, nsar, ngnss, nobs, &
                                      cny_fault, coor_fault, node_to_elem_val, node_to_elem_size, &
                                      id_dof, luni, lmat, lmat_index, lmat_val, llmat, gmat, &
                                      slip_dist, obs_points, obs_unitvec, obs_sigma, sigma2_full, &
                                      target_id_val, node_id_in_patch, xinode, etanode, uxinode, &
                                      uetanode, r1vec, r2vec, nvec, response_dist, uobs, uret, &
                                      slip_particles, slip_particles_new, nparticle_slip, max_slip, &
                                      dvec, slip_likelihood_ls, slip_prior_ls, slip_weights, slip_mean, &
                                      slip_cov, slip_likelihood_ls_new, slip_prior_ls_new, slip_assigned_num, &
                                      slip_id_start, slip_st_rand_ls, slip_metropolis_ls, gsvec, lsvec, &
                                      slip_particle_cur, slip_particle_cand, slip_st_rand)
        call mpi_barrier(mpi_comm_world, ierr)
        call mpi_gather(work_likelihood_ls, work_size, mpi_double_precision, &
                        likelihood_ls, work_size, mpi_double_precision, &
                        0, mpi_comm_world, ierr)
        if (myid == 0) then
            ! output result of stage 0(disabled)
            write (iter_char, "(i0)") iter
            filename = trim(trim(output_dir)//trim(iter_char)//".csv")
            print *, filename
            open (17, file=filename, status='replace')
            do iparticle = 1, nparticle
                do idim = 1, ndim
                    write (17, "(f12.5)", advance="no") particles(idim, iparticle)
                end do
                write (17, "(f12.5)") likelihood_ls(iparticle)
            end do
            close (17)
        end if
        iter = iter + 1
        gamma = 0d0
        do while (1.-gamma > 10d0**(-8d0))
            if (myid == 0) then
                ! find the gamma such that c.o.v of weights = 0.5
                gamma_prev = gamma
                call fault_find_next_gamma(gamma, gamma_prev, likelihood_ls, weights, nparticle)
                print *, "gamma: ", gamma
                ! normalize weights(sum of weights needs to be 1)
                call fault_normalize_weights(weights, nparticle)
                ! calculate mean and covariance of the samples
                call fault_calc_mean_particles(particles, weights, mean, nparticle, ndim)
                call fault_calc_cov_particles(particles, weights, mean, cov, &
                                              nparticle, ndim)
                call fault_resample_particles(nparticle, weights, assigned_num)
                call fault_reorder_to_send(assigned_num, particles, tmp_assigned_num, &
                                           tmp_particles, sorted_idx, nparticle, ndim, &
                                           numprocs, work_size)
                ! for(iproc=0 iproc < numprocs iproc + +) {
                do iproc = 1, numprocs
                    sum_assigned = 0
                    do iparticle = (iproc - 1)*work_size + 1, iproc*work_size
                        sum_assigned = sum_assigned + assigned_num(iparticle)
                    end do
                    sum_assigned_ls(iproc) = sum_assigned
                end do
            end if
            call mpi_bcast(gamma, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
            call mpi_bcast(cov, ndim*ndim, mpi_double_precision, 0, mpi_comm_world, ierr)
            call mpi_scatter(particles, work_size*ndim, mpi_double_precision, &
                             work_particles, work_size*ndim, mpi_double_precision, &
                             0, mpi_comm_world, ierr)
            call mpi_scatter(assigned_num, work_size, mpi_integer, work_assigned_num, &
                             work_size, mpi_integer, 0, mpi_comm_world, ierr)
!     if (myid == 0) {
!     en_time = MPI_Wtime()
!     printf("1part time: %f\n", en_time - st_time)
!     st_time = MPI_Wtime()
!     }

            id_start(1) = 1
            do iparticle = 1, work_size - 1
                id_start(iparticle + 1) = &
                    id_start(iparticle) + work_assigned_num(iparticle)
            end do
            sum_assigned = &
                id_start(work_size) + work_assigned_num(work_size) - 1
            if (sum_assigned > 50 * work_size) then
                print *, "too many samples assigned to process", myid
            end if
            call mpi_barrier(mpi_comm_world, ierr)
            st_time = omp_get_wtime()
            call work_mcmc_sampling(work_assigned_num, work_particles, &
                                    work_particles_new, &
                                    work_likelihood_ls, id_start, &
                                    particle_cur, particle_cand, &
                                    st_rand, work_size, ndim, &
                                    cov, gamma, myid, nxi, neta, nnode, ndof, &
                                    nsar, ngnss, nobs, cny_fault, coor_fault, &
                                    node_to_elem_val, node_to_elem_size, id_dof, &
                                    luni, lmat, lmat_index, lmat_val, llmat, gmat, &
                                    slip_dist, obs_points, obs_unitvec, obs_sigma, &
                                    sigma2_full, target_id_val, node_id_in_patch, &
                                    xinode, etanode, uxinode, uetanode, r1vec, r2vec, &
                                    nvec, response_dist, uobs, uret, slip_particles, &
                                    slip_particles_new, nparticle_slip, max_slip, dvec, &
                                    slip_likelihood_ls, slip_prior_ls, slip_weights, &
                                    slip_mean, slip_cov, slip_likelihood_ls_new, slip_prior_ls_new, &
                                    slip_assigned_num, slip_id_start, slip_st_rand_ls, &
                                    slip_metropolis_ls, gsvec, lsvec, slip_particle_cur, &
                                    slip_particle_cand, slip_st_rand)
            call mpi_barrier(mpi_comm_world, ierr)
            en_time = omp_get_wtime()
            if (myid == 0) then
                print *, "work_mcmc_sampling etime: ", en_time - st_time
            end if
            if (myid == 0) then
                do iproc = 2, numprocs
                    call mpi_recv(buf_likelihood(1, iproc), sum_assigned_ls(iproc), &
                                  mpi_double_precision, iproc - 1, 0, mpi_comm_world, istatus, ierr)
                    call mpi_recv(buf_particles_flat(1, iproc), &
                                  sum_assigned_ls(iproc)*ndim, mpi_double_precision, iproc - 1, 1, &
                                  mpi_comm_world, istatus, ierr)
                end do
                do iparticle = 1, sum_assigned_ls(1)
                    buf_likelihood(iparticle, 1) = work_likelihood_ls(iparticle)
                    do idim = 1, ndim
                        buf_particles_flat((iparticle - 1)*ndim + idim, 1) = &
                            work_particles_new(idim, iparticle)
                    end do
                end do
                cnt = 1
                do iproc = 1, numprocs
                    do iparticle = 1, sum_assigned_ls(iproc)
                        likelihood_ls(cnt) = buf_likelihood(iparticle, iproc)
                        do idim = 1, ndim
                            particles(idim, cnt) = &
                                buf_particles_flat((iparticle - 1)*ndim + idim, iproc)
                        end do
                        cnt = cnt + 1
                    end do
                end do
            else
                call mpi_send(work_likelihood_ls, sum_assigned, mpi_double_precision, 0, 0, &
                              mpi_comm_world, ierr)
                call mpi_send(work_particles_new, sum_assigned*ndim, mpi_double_precision, &
                              0, 1, mpi_comm_world, ierr)
            end if

            if (myid == 0) then
                ! output result of stage 0(disabled)
                write (iter_char, "(i0)") iter
                filename = trim(trim(output_dir)//trim(iter_char)//".csv")
                print *, filename
                open (17, file=filename, status='replace')
                do iparticle = 1, nparticle
                    do idim = 1, ndim
                        write (17, "(f12.5)", advance="no") particles(idim, iparticle)
                    end do
                    write (17, "(f12.5)") likelihood_ls(iparticle)
                end do
                close (17)
            end if
            iter = iter + 1
            call mpi_barrier(mpi_comm_world, ierr)
            ! if (iter == 2) then 
            !     call mpi_finalize(ierr)
            !     stop
            ! end if
!     if (myid == 0) {
!     en_time = MPI_Wtime()
!     printf("3part time: %f\n", en_time - st_time)
!     }
!     }
        end do

!     ////calculate mean and covariance of the samples
!     //if(myid == 0) {
!     //std::vector < double > weights_fin(nparticle, 1./nparticle)
!     //std::vector < double > mean =
!     //calc_mean_particles(particles_flat, weights_fin, nparticle,
!     //ndim)
!     //std::vector < double > cov_flat = calc_cov_particles(
!     //particles_flat, weights_fin, mean, nparticle, ndim)
!     //}
!     //if(myid == 0) {
!     //free(likelihood_ls)
!     //free(assigned_num)
!     //}
!     //free(work_particles_flat)
!     //free(work_init_likelihood)
!     //free(particle_cur)
!     //free(particle_cand)
!     //free(cov_flat)
!     //free(work_assigned_num)
!     return
!     }
    end subroutine fault_smc_exec
end module smc_fault
