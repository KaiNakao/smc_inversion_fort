program main
    use mpi
    use omp_lib
    use init
    use gfunc
    use smc_fault
    implicit none
    integer :: i, j
    ! number of fault planes
    integer :: nplane
    ! number of samples
    integer :: nparticle_slip, nparticle_fault
    ! number of patches
    integer :: nxi, neta
    ! dimension of  parameter
    integer :: ndim_fault, ndim_slip
    ! number of patches, nodes, dof
    integer :: npatch, nnode, ndof

    ! constrain max value for slip
    double precision :: max_slip

    ! fault geometry
    integer, allocatable :: cny_fault(:, :), &
                            node_to_elem_val(:, :), &
                            node_to_elem_size(:), &
                            id_dof(:)
    double precision, allocatable :: coor_fault(:, :)

    ! observation data
    character(len=200) :: observation_path
    double precision, allocatable :: obs_points(:, :), &
        obs_unitvec(:, :), obs_sigma(:), dvec(:)
    integer :: nsar, ngnss, nobs

    ! greens function
    double precision, allocatable :: gmat(:, :), slip_dist(:, :), &
        response_dist(:, :)
    double precision :: xinode(4), etanode(4), uxinode(4), uetanode(4), &
        r1vec(2), r2vec(2), nvec(4), uobs(3), uret(3)
    integer :: target_id_val(4), node_id_in_patch(4)

    ! fault parameter
    double precision, allocatable :: particle(:), tmp(:), slip(:, :), svec(:)

    ! laplacian
    double precision, allocatable :: luni(:, :), lmat(:, :), &
        lmat_val(:, :), llmat(:, :), ltmat_val(:, :)
    integer, allocatable :: lmat_index(:, :), ltmat_index(:, :)

    ! SMC for slip
    double precision, allocatable :: slip_particles(:, :), &
        slip_particles_new(:, :), slip_likelihood_ls(:), slip_prior_ls(:), &
        slip_weights(:), slip_mean(:), slip_cov(:, :), slip_likelihood_ls_new(:), &
        slip_prior_ls_new(:), slip_st_rand(:), slip_st_rand_ls(:, :), &
        slip_metropolis_ls(:), sigma2_full(:), alpha2_full(:), gsvec(:), lsvec(:), &
        slip_particle_cur(:), slip_particle_cand(:)
    integer, allocatable :: slip_assigned_num(:), slip_id_start(:)
    double precision :: neglog

    ! SMC for fault
    double precision, allocatable :: range(:, :)
    character(len=500) :: output_dir, command

    ! MPI
    integer :: myid, numprocs, nthreads
    integer :: ierr

    ! random seed
    integer :: seedsize
    integer, allocatable:: seed(:)

    ! measuring time
    double precision :: st_time, en_time

    double precision :: dtmp1, dtmp2

    ! set MPI
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world, myid, ierr)
    call mpi_comm_size(mpi_comm_world, numprocs, ierr)

    ! random seed
    call random_seed(size=seedsize)
    allocate (seed(seedsize))
    do i = 1, seedsize
        seed(i) = myid*100
        ! call system_clock(count=seed(i))
    end do
    call random_seed(put=seed)

    ! read settings
    if (myid == 0) then
        open (17, file="data/setting.dat", status="old")
        read (17, *) ! number of fault planes
        read (17, *) nplane
        read (17, *)
        read (17, *) ! number of samples for fault
        read (17, *) nparticle_fault
        read (17, *) ! number of samples for slip
        read (17, *) nparticle_slip
        read (17, *)
        read (17, *) ! observation path
        read (17, "(a)") observation_path
        read (17, *) ! output directory path
        read (17, "(a)") output_dir
        read (17, *)
        read (17, *) ! max slip(m)
        read (17, *) max_slip
        close (17)
        print *, "nplane: ", nplane
        print *, "nparticle_fault: ", nparticle_fault
        print *, "nparticle_slip: ", nparticle_slip
        print *, "observation_path: ", observation_path
        print *, "output_dir: ", output_dir
        print *, "max_slip: ", max_slip
    end if
    call mpi_bcast(nplane, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(nparticle_fault, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(nparticle_slip, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(max_slip, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
    if (mod(nparticle_fault, numprocs) /= 0) then
        if (myid == 0) then
            print *, "invalid number of MPI processes"
        end if
        call mpi_finalize(ierr)
        stop
    end if

    ! generate output directory
    if (myid == 0) then
        command = "mkdir -p "
        command = trim(command)//" "//output_dir
        call system(command)
    end if

    ! number of patches (per fault plane)
    nxi = 6
    neta = 6

    ! number of patches, nodes, dof
    npatch = nxi*neta*nplane
    nnode = (nxi + 1)*(neta + 1)*nplane
    ndof = (nxi - 1)*(neta - 1)*nplane

    ! dimension of fault parameter
    ndim_fault = 2 + 8*nplane
    ndim_slip = 2*ndof

    ! read observation data
    ! synthetic test
    ! observation_path = "input/observation_toy2.csv"
    ! real observation data
    ! observation_path = "input/observation_with_gnss_reduced.csv"
    if (myid == 0) then
        call read_observation1(observation_path, nobs)
        print *, "nobs: ", nobs
    end if
    call mpi_bcast(nobs, 1, mpi_integer, 0, &
                   mpi_comm_world, ierr)
    allocate (obs_points(2, nobs))
    allocate (obs_unitvec(3, nobs))
    allocate (obs_sigma(nobs))
    allocate (dvec(nobs))
    if (myid == 0) then
        call read_observation2(observation_path, nobs, &
                               obs_points, &
                               obs_unitvec, obs_sigma, &
                               dvec, nsar, ngnss)
        print *, "nsar: ", nsar
        print *, "ngnss: ", ngnss
    end if
    call mpi_bcast(nsar, 1, mpi_integer, 0, &
                   mpi_comm_world, ierr)
    call mpi_bcast(ngnss, 1, mpi_integer, 0, &
                   mpi_comm_world, ierr)
    call mpi_bcast(obs_points, 2*nobs, mpi_double_precision, &
                   0, mpi_comm_world, ierr)
    call mpi_bcast(obs_unitvec, 3*nobs, mpi_double_precision, &
                   0, mpi_comm_world, ierr)
    call mpi_bcast(obs_sigma, nobs, mpi_double_precision, &
                   0, mpi_comm_world, ierr)
    call mpi_bcast(dvec, nobs, mpi_double_precision, &
                   0, mpi_comm_world, ierr)

    ! fault geometry
    allocate (cny_fault(4, npatch))
    allocate (coor_fault(2, nnode))
    allocate (node_to_elem_val(4, nnode))
    allocate (node_to_elem_size(nnode))
    allocate (id_dof(ndof))

    ! greens function
    allocate (gmat(nobs, ndim_slip))
    allocate (slip_dist(2, nnode))
    allocate (response_dist(3, nobs))

    ! laplacian
    allocate (luni(nnode, nnode))
    allocate (lmat(2*nnode, 2*ndof))
    allocate (lmat_index(5, 2*nnode))
    allocate (lmat_val(5, 2*nnode))
    allocate (ltmat_index(5, 2*ndof))
    allocate (ltmat_val(5, 2*ndof))
    allocate (llmat(2*ndof, 2*ndof))

    ! SMC for slip
    allocate (slip_particles(ndim_slip, nparticle_slip))
    allocate (slip_particles_new(ndim_slip, nparticle_slip))
    allocate (slip_likelihood_ls(nparticle_slip))
    allocate (slip_prior_ls(nparticle_slip))
    allocate (slip_weights(nparticle_slip))
    allocate (slip_mean(ndim_slip))
    allocate (slip_cov(ndim_slip, ndim_slip))
    allocate (slip_likelihood_ls_new(nparticle_slip))
    allocate (slip_prior_ls_new(nparticle_slip))
    allocate (slip_st_rand(ndim_slip))
    allocate (slip_assigned_num(nparticle_slip))
    allocate (slip_id_start(nparticle_slip))
    allocate (slip_st_rand_ls(ndim_slip, nparticle_slip))
    allocate (slip_metropolis_ls(nparticle_slip))
    allocate (sigma2_full(nobs))
    allocate (alpha2_full(2*nnode))
    allocate (gsvec(nobs))
    allocate (lsvec(2*nnode))
    allocate (slip_particle_cur(ndim_slip))
    allocate (slip_particle_cand(ndim_slip))

    ! ! calculate diplacement for fixed fault and slip distribution
    ! allocate (slip(2, nnode))
    ! allocate (svec(2*ndof))

    ! allocate (particle(ndim_fault))
    ! open (10, file="mean_fault.dat", status='old')
    ! ! open (10, file="/home/nakao/smc_inversion_fort/input/noto_synthetic/theta.dat", &
    ! !       status='old')
    ! ! open (10, file="data/theta.dat")
    ! do i = 1, ndim_fault
    !     read (10, *) particle(i)
    ! end do
    ! close (10)

    ! print *, particle
    ! call discretize_fault(particle, nplane, nxi, neta, cny_fault, coor_fault, &
    !                       node_to_elem_val, node_to_elem_size, id_dof)
    ! open (10, file="mean_slip.dat", &
    !       status='old')
    ! do i = 1, nnode
    !     read (10, *) slip(1, i), slip(2, i)
    ! end do
    ! do i = 1, ndof
    !     svec(2*i - 1) = slip(1, id_dof(i))
    !     svec(2*i) = slip(2, id_dof(i))
    ! end do
    ! close (10)
    ! print *, "svec:", svec

    ! call calc_greens_func(particle, nplane, nxi, neta, gmat, slip_dist, cny_fault, coor_fault, obs_points, &
    !                       obs_unitvec, node_to_elem_val, node_to_elem_size, &
    !                       id_dof, nsar, ngnss, nobs, nnode, ndof, target_id_val, &
    !                       node_id_in_patch, xinode, etanode, uxinode, uetanode, &
    !                       r1vec, r2vec, nvec, response_dist, uobs, uret)

    ! call dgemv('n', nobs, 2*ndof, 1d0, gmat, &
    !            nobs, svec, 1, 0d0, gsvec, 1)
    ! open (10, file="dvec_est.dat", &
    !       status='replace')
    ! ! open (10, file="/home/nakao/smc_inversion_fort/input/noto_synthetic/dvec_exact.dat", &
    ! !       status='replace')
    ! do i = 1, nobs
    !     write (10, *), gsvec(i)
    ! end do
    ! close (10)
    ! print *, "gsvec, dvec"

    ! dtmp1 = 0d0
    ! dtmp2 = 0d0
    ! do i = 1, nobs
    !     dtmp1 = dtmp1 + (gsvec(i) - dvec(i))**2/obs_sigma(i)**2
    !     dtmp2 = dtmp2 + dvec(i)**2/obs_sigma(i)**2
    !     write (*, *), gsvec(i), dvec(i)
    ! end do
    ! print *, "variance reduction: ", 1d0 - dtmp1/dtmp2

    ! ! calculate likelihood for given fault
    ! allocate (particle(ndim_fault))
    ! allocate (tmp(ndim_fault + 1))
    ! do i = 1, ndim_fault
    !     particle(i) = 0d0
    ! end do
    ! do i = 1, ndim_fault + 1
    !     tmp(i) = 0d0
    ! end do
    ! open (10, file="64.csv", status='old')
    ! ! open (10, file="output_obs_1e4_1e3_3plane_gaussian/57.csv", status='old')
    ! do i = 1, nparticle_fault
    !     ! read (10, *) tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), &
    !     !     tmp(6), tmp(7), tmp(8), tmp(9), tmp(10), tmp(11)
    !     read (10, *) tmp
    !     do j = 1, ndim_fault
    !         particle(j) = particle(j) + tmp(j)
    !     end do
    ! end do
    ! close (10)
    ! tmp = particle
    ! do j = 1, ndim_fault
    !     particle(j) = particle(j)/nparticle_fault
    ! end do
    ! open (10, file="mean_fault.dat", status="replace")
    ! do i = 1, ndim_fault
    !     write (10, *) particle(i)
    ! end do
    ! close (10)

    ! ! ! ! ! ! open (10, file="/hoe/nakao/smc_inversion_fort/input/noto_synthetic/theta.dat", &
    ! ! ! ! ! !       status="old")
    ! ! ! open (10, file="data/theta.dat", status="old")
    ! ! open (10, file="mean_fault.dat", status="old")
    ! ! do i = 1, ndim_fault
    ! !     read (10, *) particle(i)
    ! ! end do
    ! ! close (10)
    ! ! print *, particle

    ! st_time = omp_get_wtime()
    ! print *, "start"
    ! neglog = fault_calc_likelihood( &
    !          particle, nplane, nxi, neta, nnode, ndof, nsar, ngnss, nobs, cny_fault, &
    !          coor_fault, node_to_elem_val, node_to_elem_size, id_dof, luni, lmat, &
    !          lmat_index, lmat_val, ltmat_index, ltmat_val, llmat, gmat, slip_dist, obs_points, &
    !          obs_unitvec, obs_sigma, sigma2_full, alpha2_full, target_id_val, node_id_in_patch, &
    !          xinode, etanode, uxinode, uetanode, r1vec, r2vec, nvec, response_dist, &
    !          uobs, uret, slip_particles, slip_particles_new, &
    !          nparticle_slip, max_slip, dvec, slip_likelihood_ls, slip_prior_ls, &
    !          slip_weights, slip_mean, slip_cov, slip_likelihood_ls_new, &
    !          slip_prior_ls_new, slip_assigned_num, slip_id_start, slip_st_rand_ls, &
    !          slip_metropolis_ls, gsvec, lsvec, slip_particle_cur, &
    !          slip_particle_cand, slip_st_rand, 1, "output/slip_from_mean_fault.dat")
    ! !  slip_particle_cand, slip_st_rand, 0, "output/slip_from_mean_fault.dat")
    ! en_time = omp_get_wtime()
    ! print *, "etime: ", en_time - st_time
    ! print *, "neglog: ", neglog

    ! call calc_slip_map("64.csv")
    if (myid == 0) then
        call point_cloud("64.csv", &
                         "output/mapslip.dat")
    end if

    ! ! smc for fault
    ! allocate (range(2, ndim_fault))
    ! ! range of uniform prior distribution P(theta)
    ! ! ! synthetic test
    ! ! ! range(:, :) = reshape((/-5., 15., -15., 15., -39., -10., -20., 20., 50., 90., &
    ! ! !                         -2., 2., -2., 2., -10., 2., 1., 50., 1., 50./), &
    ! ! !                       (/2, ndim_fault/))
    ! ! ! real observation data
    ! ! range(:, :) = reshape((/-10, 10, -30, 0, -30, -1, -20, 20, 50, 90, &
    ! !                         -2, 2, -2, 2, -10, 2, 1, 50, 1, 50/), &
    ! !                       (/2, ndim_fault/))

    ! ! read prior range of theta
    ! if (myid == 0) then
    !     open (17, file="data/range.dat", status="old")
    !     do i = 1, ndim_fault
    !         read (17, "(a)") ! param
    !         read (17, *) range(1, i), range(2, i)
    !     end do
    !     close (17)
    !     print *, "prior range"
    !     do i = 1, nplane
    !         print *, "fault plane: ", i
    !         print *, "xf ", range(1, 8*i - 7), range(2, 8*i - 7)
    !         print *, "yf ", range(1, 8*i - 6), range(2, 8*i - 6)
    !         print *, "zf ", range(1, 8*i - 5), range(2, 8*i - 5)
    !         print *, "strike ", range(1, 8*i - 4), range(2, 8*i - 4)
    !         print *, "dip ", range(1, 8*i - 3), range(2, 8*i - 3)
    !         print *, "lxi ", range(1, 8*i - 2), range(2, 8*i - 2)
    !         print *, "leta ", range(1, 8*i - 1), range(2, 8*i - 1)
    !         print *, "log_alpha2 ", range(1, 8*i - 0), range(2, 8*i - 0)
    !     end do
    !     print *, "log_sigma_sar2 ", range(1, 8*nplane + 1), range(2, 8*nplane + 1)
    !     print *, "log_sigma_gnss2 ", range(1, 8*nplane + 2), range(2, 8*nplane + 2)
    ! end if
    ! call mpi_bcast(range, 2*ndim_fault, mpi_double_precision, 0, &
    !                mpi_comm_world, ierr)

    ! call fault_smc_exec( &
    !     output_dir, range, nplane, nparticle_fault, ndim_fault, &
    !     myid, numprocs, nxi, neta, nnode, ndof, nsar, ngnss, nobs, cny_fault, &
    !     coor_fault, node_to_elem_val, node_to_elem_size, id_dof, luni, lmat, &
    !     lmat_index, lmat_val, ltmat_index, ltmat_val, llmat, gmat, slip_dist, obs_points, &
    !     obs_unitvec, obs_sigma, sigma2_full, alpha2_full, target_id_val, node_id_in_patch, &
    !     xinode, etanode, uxinode, uetanode, r1vec, r2vec, nvec, &
    !     response_dist, uobs, uret, slip_particles, &
    !     slip_particles_new, nparticle_slip, max_slip, dvec, gsvec, &
    !     lsvec, slip_likelihood_ls, slip_prior_ls, slip_weights, slip_mean, &
    !     slip_cov, slip_likelihood_ls_new, slip_prior_ls_new, &
    !     slip_st_rand, slip_particle_cur, slip_particle_cand, &
    !     slip_assigned_num, slip_id_start, slip_st_rand_ls, slip_metropolis_ls)
    call mpi_finalize(ierr)

contains
    subroutine calc_slip_map(fault_filepath)
        implicit none
        character(*), intent(in) :: fault_filepath
        character(len=100) :: slip_filepath
        double precision, allocatable :: fault_particles(:, :), &
            work_fault_particles(:, :), mean_slip(:, :), &
            work_mean_slip(:, :), particle_slip(:)
        integer :: work_size, k, inode, idirection, iparticle, idim, idof
        integer :: stride

        stride = 100

        slip_filepath = "output/mapslip.dat"
        allocate (slip(2, nnode))
        allocate (particle_slip(ndim_slip))
        if (myid == 0) then
            allocate (fault_particles(ndim_fault, nparticle_fault))
            allocate (mean_slip(ndim_slip, nparticle_fault))
            allocate (tmp(ndim_fault + 1))
            open (10, file=fault_filepath, status='old')
            do i = 1, nparticle_fault
                do k = 1, stride
                    read (10, *) tmp
                end do
                do j = 1, ndim_fault
                    fault_particles(j, i) = tmp(j)
                end do
            end do
            close (10)
        else
            allocate (fault_particles(1, 1))
            allocate (mean_slip(1, 1))
        end if

        work_size = nparticle_fault/numprocs
        allocate (work_fault_particles(ndim_fault, work_size))
        allocate (work_mean_slip(ndim_slip, work_size))
        call mpi_scatter(fault_particles, work_size*ndim_fault, &
                         mpi_double_precision, work_fault_particles, work_size*ndim_fault, &
                         mpi_double_precision, 0, mpi_comm_world, ierr)

        allocate (particle(ndim_fault))
        do i = 1, work_size
            do j = 1, ndim_fault
                particle(j) = work_fault_particles(j, i)
            end do

            neglog = fault_calc_likelihood(particle, nplane, nxi, neta, nnode, ndof, &
                                           nsar, ngnss, nobs, cny_fault, coor_fault, node_to_elem_val, &
                                           node_to_elem_size, id_dof, luni, lmat, lmat_index, lmat_val, &
                                           ltmat_index, ltmat_val, llmat, gmat, slip_dist, obs_points, &
                                           obs_unitvec, obs_sigma, sigma2_full, alpha2_full, target_id_val, &
                                           node_id_in_patch, xinode, etanode, uxinode, uetanode, &
                                           r1vec, r2vec, nvec, response_dist, uobs, uret, slip_particles, &
                                           slip_particles_new, nparticle_slip, max_slip, dvec, slip_likelihood_ls, &
                                           slip_prior_ls, slip_weights, slip_mean, slip_cov, slip_likelihood_ls_new, &
                                           slip_prior_ls_new, slip_assigned_num, slip_id_start, slip_st_rand_ls, &
                                           slip_metropolis_ls, gsvec, lsvec, slip_particle_cur, slip_particle_cand, &
                                           slip_st_rand, 0, "")

            print *, "ID: ", work_size*myid + i
            do j = 1, ndim_slip
                work_mean_slip(j, i) = 0d0
            end do

            do k = 1, nparticle_slip
                do j = 1, ndim_slip
                    work_mean_slip(j, i) = &
                        work_mean_slip(j, i) + slip_particles(j, k)
                end do
            end do
            do j = 1, ndim_slip
                work_mean_slip(j, i) = work_mean_slip(j, i)/nparticle_slip
            end do
        end do

        call mpi_gather(work_mean_slip, ndim_slip*work_size, &
                        mpi_double_precision, mean_slip, ndim_slip*work_size, &
                        mpi_double_precision, 0, mpi_comm_world, ierr)

        if (myid == 0) then
            open (17, file=slip_filepath, status='replace')
            do iparticle = 1, nparticle_fault
                do inode = 1, nnode
                    do idirection = 1, 2
                        slip(idirection, inode) = 0d0
                    end do
                end do
                do idim = 1, ndim_slip
                    particle_slip(idim) = mean_slip(idim, iparticle)
                end do
                do idof = 1, ndof
                    inode = id_dof(idof)
                    slip(1, inode) = particle_slip(2*(idof - 1) + 1)
                    slip(2, inode) = particle_slip(2*(idof - 1) + 2)
                end do
                do inode = 1, nnode
                    do idirection = 1, 2
                        write (17, "(f12.5)", advance="no") slip(idirection, inode)
                    end do
                end do
                write (17, *)
            end do
            close (17)
        end if
    end subroutine calc_slip_map

    subroutine point_cloud(fault_filepath, slip_filepath)
        implicit none
        character(*), intent(in) :: fault_filepath, slip_filepath
        character(len=200) :: points_filepath
        double precision, allocatable :: fault_particles(:, :), mean_slip(:, :), &
            points(:, :)
        double precision :: r1, r2, xi, eta, uxi, ueta, x, y, z, &
            xf, yf, zf, strike, dip, lxi, leta, dip_rad, strike_rad, &
            theta(ndim_fault)

        integer :: ie, iparticle, idim, inode, ir1, ir2, iplane, jnode, idpoint, cnt, stride, k
        double precision :: pi = 4d0*atan(1d0)

        integer :: nr = 4, ir
        double precision, allocatable :: r1vec_cloud(:), r2vec_cloud(:)

        allocate(r1vec_cloud(nr))
        allocate(r2vec_cloud(nr))
        do ir = 1, nr
            r1vec_cloud(ir) = -1d0 + 2d0 / dble(nr)
            r2vec_cloud(ir) = -1d0 + 2d0 / dble(nr)
        end do

        stride = 100
        points_filepath = "output/points.dat"
        allocate (fault_particles(ndim_fault, nparticle_fault))
        ! x, y, z, uxi, ueta, unorm
        allocate (points(6, nparticle_fault*npatch*nr**2))
        allocate (mean_slip(2*nnode, nparticle_fault))
        allocate (tmp(ndim_fault + 1))
        allocate (particle(ndim_fault))
        allocate (svec(ndim_slip))
        open (10, file=fault_filepath, status='old')
        open (11, file=slip_filepath, status='old')
        do i = 1, nparticle_fault
            do k = 1, stride
                read (10, *) tmp
            end do
            read (11, *) mean_slip(:, i)
            do j = 1, ndim_fault
                fault_particles(j, i) = tmp(j)
            end do
        end do
        close (10)
        close (11)

        cnt = 1
        do iparticle = 1, nparticle_fault
            print *, "iparticle: ", iparticle
            slip_dist = reshape(mean_slip(:, iparticle), (/2, nnode/))
            do idim = 1, ndim_fault
                particle(idim) = fault_particles(idim, iparticle)
            end do

            call discretize_fault(particle, nplane, nxi, neta, cny_fault, coor_fault, &
                                  node_to_elem_val, node_to_elem_size, id_dof)

            ! loop for patchs
            do iplane = 1, nplane
                xf = particle(8*(iplane-1)+1)
                yf = particle(8*(iplane-1)+2)
                zf = particle(8*(iplane-1)+3)
                strike = particle(8*(iplane-1)+4)
                strike_rad = strike*pi/180d0
                dip = particle(8*(iplane-1)+5)
                dip_rad = dip*pi/180d0
                lxi = particle(8*(iplane-1)+6)
                leta = particle(8*(iplane-1)+7)
                do ie = 1 + nxi*neta*(iplane - 1), nxi*neta*iplane
                    do inode = 1, 4
                        node_id_in_patch(inode) = cny_fault(inode, ie)
                    end do

                    do inode = 1, 4
                        jnode = node_id_in_patch(inode)
                        xinode(inode) = coor_fault(1, jnode)
                        etanode(inode) = coor_fault(2, jnode)
                        uxinode(inode) = slip_dist(1, jnode)
                        uetanode(inode) = slip_dist(2, jnode)
                    end do

                    !   loop for point sources
                    do ir1 = 1, nr
                        r1 = r1vec_cloud(ir1)
                        do ir2 = 1, nr
                            r2 = r2vec_cloud(ir2)
                            !   interpolate(xi, eta, uxi, ueta) by shape functions
                            nvec(1) = (1d0 - r1)*(1d0 - r2)/4d0
                            nvec(2) = (1d0 + r1)*(1d0 - r2)/4d0
                            nvec(3) = (1d0 + r1)*(1d0 + r2)/4d0
                            nvec(4) = (1d0 - r1)*(1d0 + r2)/4d0
                            xi = 0d0
                            eta = 0d0
                            uxi = 0d0
                            ueta = 0d0
                            do inode = 1, 4
                                xi = xi + nvec(inode)*xinode(inode)
                                eta = eta + nvec(inode)*etanode(inode)
                                uxi = uxi + nvec(inode)*uxinode(inode)
                                ueta = ueta + nvec(inode)*uetanode(inode)
                            end do

                            !   location of the point source specified by
                            !   global coordinates(x, y, z)
                            x = xf - eta*cos(dip_rad)*cos(strike_rad) + &
                                xi*sin(strike_rad)
                            y = yf + eta*cos(dip_rad)*sin(strike_rad) + &
                                xi*cos(strike_rad)
                            z = zf + eta*sin(dip_rad)
                            idpoint = (cnt - 1)*npatch*4 + (ie - 1)*4 + (ir1 - 1)*2 + ir2
                            points(1, idpoint) = x
                            points(2, idpoint) = y
                            points(3, idpoint) = z
                            points(4, idpoint) = uxi
                            points(5, idpoint) = ueta
                            points(6, idpoint) = sqrt(uxi**2 + ueta**2)
                        end do
                    end do
                end do
            end do
            cnt = cnt + 1
        end do

        open (10, file=points_filepath, status='replace')
        do i = 1, nparticle_fault*npatch*nr**2
            do j = 1, 6
                write (10, "(f12.5)", advance="no") points(j, i)
            end do
            write (10, *)
        end do
        close (10)
    end subroutine point_cloud

end program main
