program main
    use type_mat
    use mpi
    use omp_lib
    use init
    use gfunc
    use smc_fault
    implicit none
    integer :: i, j
    ! number of fault planes
    integer :: iplane, nplane
    ! number of samples
    integer :: nparticle_slip, nparticle_fault
    ! number of patches
    integer :: nxi, neta
    integer, allocatable :: nxi_ls(:), neta_ls(:)
    ! dimension of  parameter
    integer :: ndim_fault, ndim_slip
    ! number of patches, nodes, dof (per fault plane)
    integer, allocatable :: npatch_ls(:), nnode_ls(:), ndof_ls(:), ndof_index(:)
    integer :: npatch_total, nnode_total, ndof_total

    ! constrain max value for slip
    double precision :: max_slip

    ! fault geometry
    integer, allocatable :: cny_fault(:, :), &
                            node_to_elem_val(:, :), &
                            node_to_elem_size(:), &
                            id_dof(:)
    double precision, allocatable :: coor_fault(:, :)
    double precision :: xmin, xmax, zmin
    logical :: fix_xbend
    double precision, allocatable :: xbend(:)

    ! observation data
    character(len=200) :: observation_path
    double precision, allocatable :: obs_points(:, :), &
        obs_unitvec(:, :), obs_sigma(:), dvec(:)
    integer :: npath
    integer, allocatable :: nsar_ls(:), nsar_index(:)
    integer :: nsar_total, ngnss, nobs

    ! greens function
    double precision, allocatable :: gmat(:, :), slip_dist(:, :), &
        response_dist(:, :)
    type(mat), allocatable :: gmat_arr(:)
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
    type(mat), allocatable :: sigma_sar_mat(:)

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

    ! Generalized SMC
    logical :: gsmc
    double precision, allocatable :: qmean(:), qcov(:, :), qcov_inv(:, :)
    double precision :: log_qcov_det

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
        allocate (nxi_ls(nplane))
        allocate (neta_ls(nplane))
        allocate (npatch_ls(nplane))
        allocate (nnode_ls(nplane))
        allocate (ndof_ls(nplane))
        allocate (ndof_index(nplane + 1))
        allocate (xbend(nplane - 1))
        read (17, *) ! fix-xbend
        read (17, *) fix_xbend
        read (17, *) ! xbend
        read (17, *) xbend
        read (17, *)
        read (17, *) ! number of samples for fault
        read (17, *) nparticle_fault
        read (17, *) ! number of samples for slip
        read (17, *) nparticle_slip
        read (17, *)
        read (17, *) ! observation path
        read (17, "(a)") observation_path
        read (17, *) ! number of SAR paths
        read (17, *) npath
        allocate (nsar_ls(npath))
        allocate (nsar_index(npath + 1))
        allocate (sigma_sar_mat(npath))
        allocate (gmat_arr(npath))
        read (17, *) ! output directory path
        read (17, "(a)") output_dir
        read (17, *)
        read (17, *) ! max slip(m)
        read (17, *) max_slip
        read (17, *)
        read (17, *) ! nxi
        read (17, *) nxi_ls
        read (17, *) ! neta
        read (17, *) neta_ls
        read (17, *)
        read (17, *) ! xmin
        read (17, *) xmin
        read (17, *) ! xmax
        read (17, *) xmax
        read (17, *) ! zmin
        read (17, *) zmin
        read (17, *)
        read (17, *) ! GSMC
        read (17, *) gsmc
        close (17)
        print *, "nplane: ", nplane
        print *, "nparticle_fault: ", nparticle_fault
        print *, "nparticle_slip: ", nparticle_slip
        print *, "observation_path: ", observation_path
        print *, "output_dir: ", output_dir
        print *, "max_slip: ", max_slip
        print *, "nxi_ls: ", nxi_ls
        print *, "neta_ls: ", neta_ls
    end if
    call mpi_bcast(nplane, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(fix_xbend, 1, mpi_logical, 0, mpi_comm_world, ierr)
    call mpi_bcast(npath, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(nparticle_fault, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(nparticle_slip, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(max_slip, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
    call mpi_bcast(xmin, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
    call mpi_bcast(xmax, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
    call mpi_bcast(zmin, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
    call mpi_bcast(gsmc, 1, mpi_logical, 0, mpi_comm_world, ierr)
    if (myid /= 0) then
        allocate (nxi_ls(nplane))
        allocate (neta_ls(nplane))
        allocate (nsar_ls(npath))
        allocate (nsar_index(npath + 1))
        allocate (sigma_sar_mat(npath))
        allocate (gmat_arr(npath))
        allocate (npatch_ls(nplane))
        allocate (nnode_ls(nplane))
        allocate (ndof_ls(nplane))
        allocate (ndof_index(nplane + 1))
        allocate (xbend(nplane - 1))
    end if
    call mpi_bcast(nxi_ls, nplane, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(neta_ls, nplane, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(xbend, nplane - 1, mpi_double_precision, 0, mpi_comm_world, ierr)

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

    ! number of patches, nodes, dof
    npatch_ls = 0
    npatch_total = 0
    nnode_ls = 0
    nnode_total = 0
    ndof_ls = 0
    ndof_total = 0
    ndof_index(1) = 1
    do iplane = 1, nplane
        nxi = nxi_ls(iplane)
        neta = neta_ls(iplane)
        npatch_ls(iplane) = nxi*neta
        nnode_ls(iplane) = (nxi + 1)*(neta + 1)
        ! if (iplane == 1 .or. iplane == nplane) then
        !     ndof_ls(iplane) = nnode_ls(iplane) - (nxi + neta + 1)
        ! else
        !     ndof_ls(iplane) = nnode_ls(iplane) - (nxi + 1)
        ! end if
        ndof_ls(iplane) = nnode_ls(iplane) - (nxi + 1)
        if (iplane == 1) then
            ndof_ls(iplane) = ndof_ls(iplane) - neta
        end if
        if (iplane == nplane) then
            ndof_ls(iplane) = ndof_ls(iplane) - neta
        end if
        npatch_total = npatch_total + npatch_ls(iplane)
        nnode_total = nnode_total + nnode_ls(iplane)
        ndof_total = ndof_total + ndof_ls(iplane)
        ndof_index(iplane + 1) = ndof_index(iplane) + ndof_ls(iplane)
    end do

    ! dimension of fault parameter
    if (fix_xbend) then
        ndim_fault = 2*nplane + 4
    else
        ndim_fault = 3*nplane + 3
    end if
    ndim_slip = 2*ndof_total
    if (myid == 0) then
        print *, "ndim_fault: ", ndim_fault
        print *, "ndim_slip: ", ndim_slip
        print *, "fix_xbend: ", fix_xbend
        print *, "xbend: ", xbend
    end if

    ! read observation data
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
                               dvec, nsar_ls, nsar_index, nsar_total, ngnss, npath)
        print *, "nsar: ", nsar_ls
        print *, "ngnss: ", ngnss
    end if
    call mpi_bcast(nsar_total, 1, mpi_integer, 0, &
                   mpi_comm_world, ierr)
    call mpi_bcast(nsar_ls, npath, mpi_integer, 0, &
                   mpi_comm_world, ierr)
    call mpi_bcast(nsar_index, npath + 1, mpi_integer, 0, &
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

    do i = 1, npath
        sigma_sar_mat(i)%ncol = nsar_ls(i)
        sigma_sar_mat(i)%nrow = nsar_ls(i)
        allocate (sigma_sar_mat(i)%body(sigma_sar_mat(i)%nrow, sigma_sar_mat(i)%ncol))

        gmat_arr(i)%ncol = ndim_slip
        gmat_arr(i)%nrow = nsar_ls(i)
        allocate (gmat_arr(i)%body(gmat_arr(i)%nrow, gmat_arr(i)%ncol))
    end do

    if (myid == 0) then
        call calc_sigma_sar_mat(npath, nsar_ls, nsar_index, &
                                obs_points, obs_sigma, sigma_sar_mat)
    end if
    do i = 1, npath
        call mpi_bcast(sigma_sar_mat(i)%body, nsar_ls(i)**2, mpi_double_precision, &
                       0, mpi_comm_world, ierr)
    end do

    ! read parameteres for distribution q
    allocate (qmean(ndim_fault))
    allocate (qcov(ndim_fault, ndim_fault))
    allocate (qcov_inv(ndim_fault, ndim_fault))
    if (myid == 0) then
        if (gsmc) then
            ! qmean
            open (10, file="data/qmean.dat", status='old')
            read (10, *) qmean
            close (10)
            ! qcov
            open (10, file="data/qcov.dat", status='old')
            read (10, *) qcov
            close (10)
            ! qcov_inv
            open (10, file="data/qcov_inv.dat", status='old')
            read (10, *) qcov_inv
            close (10)
            ! log_qcov_det
            open (10, file="data/log_qcov_det.dat", status='old')
            read (10, *) log_qcov_det
            close (10)
            ! cholesky
            call dpotrf('L', ndim_fault, qcov, ndim_fault, ierr)
        else
            qmean = 0d0
            qcov = 0d0
            qmean = 0d0
            log_qcov_det = 0d0
        end if
    end if
    call mpi_bcast(qmean, ndim_fault, mpi_double_precision, 0, &
                   mpi_comm_world, ierr)
    call mpi_bcast(qcov, ndim_fault**2, mpi_double_precision, 0, &
                   mpi_comm_world, ierr)
    call mpi_bcast(qcov_inv, ndim_fault**2, mpi_double_precision, 0, &
                   mpi_comm_world, ierr)
    call mpi_bcast(log_qcov_det, 1, mpi_double_precision, 0, &
                   mpi_comm_world, ierr)

    ! fault geometry
    allocate (cny_fault(4, npatch_total))
    allocate (coor_fault(2, nnode_total))
    allocate (node_to_elem_val(4, nnode_total))
    allocate (node_to_elem_size(nnode_total))
    allocate (id_dof(ndof_total))

! greens function
    allocate (gmat(nobs, ndim_slip))
    allocate (slip_dist(2, nnode_total))
    allocate (response_dist(3, nobs))

! laplacian
    allocate (luni(nnode_total, nnode_total))
    allocate (lmat(2*nnode_total, 2*ndof_total))
    allocate (lmat_index(9, 2*nnode_total))
    allocate (lmat_val(9, 2*nnode_total))
    allocate (ltmat_index(9, 2*ndof_total))
    allocate (ltmat_val(9, 2*ndof_total))
    allocate (llmat(2*ndof_total, 2*ndof_total))

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
    allocate (alpha2_full(2*nnode_total))
    allocate (gsvec(nobs))
    allocate (lsvec(2*nnode_total))
    allocate (slip_particle_cur(ndim_slip))
    allocate (slip_particle_cand(ndim_slip))

! ! calculate diplacement for fixed fault and slip distribution
! allocate (slip(2, nnode_total))
! allocate (svec(2*ndof_total))

! allocate (particle(ndim_fault))
! open (10, file="data/theta.dat", status='old')
! do i = 1, ndim_fault
!     read (10, *) particle(i)
! end do
! close (10)

! print *, particle
! call discretize_fault(particle, nplane, nxi_ls, neta_ls, cny_fault, coor_fault, &
!                       node_to_elem_val, node_to_elem_size, id_dof, xmin, xmax, zmin)
! open (10, file="mean_slip.dat", status='old')
! ! open (10, file="/home/nakao/smc_inversion_fort/input/noto_synthetic2/slip_toy.csv", &
! !       status='old')
! do i = 1, nnode_total
!     read (10, *) slip(1, i), slip(2, i)
! end do
! do i = 1, ndof_total
!     svec(2*i - 1) = slip(1, id_dof(i))
!     svec(2*i) = slip(2, id_dof(i))
! end do
! close (10)
! print *, "svec:", svec

! print *, "calc_greens_func start"
! call calc_greens_func(particle, nplane, nxi_ls, neta_ls, gmat, slip_dist, cny_fault, coor_fault, obs_points, &
!                       obs_unitvec, node_to_elem_val, node_to_elem_size, &
!                       id_dof, ngnss, nobs, nnode_total, ndof_total, ndof_index, target_id_val, &
!                       node_id_in_patch, xinode, etanode, uxinode, uetanode, &
!                       r1vec, r2vec, nvec, response_dist, uobs, uret, nsar_total, npath, nsar_index, gmat_arr, &
!                       xmin, xmax, zmin)
! print *, "dgemv start"
! call dgemv('n', nobs, 2*ndof_total, 1d0, gmat, &
!            nobs, svec, 1, 0d0, gsvec, 1)
! open (10, file="dvec_est.dat", status='replace')
! ! open (10, file="/home/nakao/smc_inversion_fort/input/noto_synthetic2/dvec_exact.dat", &
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

    ! calculate likelihood for given fault
    allocate (particle(ndim_fault))
    allocate (tmp(ndim_fault + 1))
    do i = 1, ndim_fault
        particle(i) = 0d0
    end do
    do i = 1, ndim_fault + 1
        tmp(i) = 0d0
    end do

    open (10, file="peak_fault.dat", status="old")
    do i = 1, ndim_fault
        read (10, *) particle(i)
    end do
    close (10)
    print *, particle

    st_time = omp_get_wtime()
    print *, "start"
    neglog = fault_calc_likelihood( &
             particle, nplane, nxi_ls, neta_ls, nnode_total, ndof_total, ndof_index, ngnss, nobs, cny_fault, &
             coor_fault, node_to_elem_val, node_to_elem_size, id_dof, luni, lmat, &
             lmat_index, lmat_val, ltmat_index, ltmat_val, llmat, gmat, slip_dist, obs_points, &
             obs_unitvec, obs_sigma, sigma2_full, alpha2_full, target_id_val, node_id_in_patch, &
             xinode, etanode, uxinode, uetanode, r1vec, r2vec, nvec, response_dist, &
             uobs, uret, slip_particles, slip_particles_new, &
             nparticle_slip, max_slip, dvec, slip_likelihood_ls, slip_prior_ls, &
             slip_weights, slip_mean, slip_cov, slip_likelihood_ls_new, &
             slip_prior_ls_new, slip_assigned_num, slip_id_start, slip_st_rand_ls, &
             slip_metropolis_ls, gsvec, lsvec, slip_particle_cur, &
             slip_particle_cand, slip_st_rand, 1, "output/slip_from_mean_fault.dat", &
             npath, nsar_index, nsar_total, sigma_sar_mat, gmat_arr, xmin, xmax, zmin, fix_xbend, xbend)
    en_time = omp_get_wtime()
    print *, "etime: ", en_time - st_time
    print *, "neglog: ", neglog
    print *, "Ndim slip: ", ndim_slip
    print *, "ABIC: ", 2d0*neglog + 2d0*ndim_fault

! call calc_slip_map("from_fugaku/shuf/38400_cov_pos.csv")
! if (myid == 0) then
!     call point_cloud("from_fugaku/shuf/192000_cov_pos.csv", &
!                      "output/shuf/mapslip_192000.dat")
! end if

    ! ! smc for fault
    ! allocate (range(2, ndim_fault))
    ! ! read prior range of theta
    ! if (myid == 0) then
    !     open (17, file="data/range.dat", status="old")
    !     do i = 1, ndim_fault
    !         read (17, "(a)") ! param
    !         read (17, *) range(1, i), range(2, i)
    !     end do
    !     close (17)

    !     print *, "prior range"
    !     print *, "ymin ", range(1, 1), range(2, 1)
    !     print *, "ymax ", range(1, 2), range(2, 2)
    !     if (fix_xbend) then
    !         if (nplane > 1) then
    !             do i = 1, nplane - 1
    !                 print *, "y", i, " ", range(1, i + 2), range(2, i + 2)
    !             end do
    !         end if
    !         do i = 1, nplane
    !             print *, "dip", i, " ", range(1, 2*nplane - 2 + i), range(2, 2*nplane - 2 + i)
    !         end do
    !         print *, "log_alpha2 ", range(1, 2*nplane + 2), range(2, 2*nplane + 2)
    !         print *, "log_sigma_sar2 ", range(1, 2*nplane + 3), range(2, 2*nplane + 3)
    !         print *, "log_sigma_gnss2 ", range(1, 2*nplane + 4), range(2, 2*nplane + 4)
    !     else
    !         if (nplane > 1) then
    !             do i = 1, nplane - 1
    !                 print *, "x", i, " ", range(1, 2*i + 1), range(2, 2*i + 1)
    !                 print *, "y", i, " ", range(1, 2*i + 2), range(2, 2*i + 2)
    !             end do
    !         end if
    !         do i = 1, nplane
    !             print *, "dip", i, " ", range(1, 2*nplane + i), range(2, 2*nplane + i)
    !         end do
    !         print *, "log_alpha2 ", range(1, 3*nplane + 1), range(2, 3*nplane + 1)
    !         print *, "log_sigma_sar2 ", range(1, 3*nplane + 2), range(2, 3*nplane + 2)
    !         print *, "log_sigma_gnss2 ", range(1, 3*nplane + 3), range(2, 3*nplane + 3)
    !     end if
    ! end if
    ! call mpi_bcast(range, 2*ndim_fault, mpi_double_precision, 0, &
    !                mpi_comm_world, ierr)

    ! call fault_smc_exec( &
    !     output_dir, range, nplane, nparticle_fault, ndim_fault, &
    !     myid, numprocs, nxi_ls, neta_ls, nnode_total, ndof_total, ndof_index, ngnss, nobs, cny_fault, &
    !     coor_fault, node_to_elem_val, node_to_elem_size, id_dof, luni, lmat, &
    !     lmat_index, lmat_val, ltmat_index, ltmat_val, llmat, gmat, slip_dist, obs_points, &
    !     obs_unitvec, obs_sigma, sigma2_full, alpha2_full, target_id_val, node_id_in_patch, &
    !     xinode, etanode, uxinode, uetanode, r1vec, r2vec, nvec, &
    !     response_dist, uobs, uret, slip_particles, &
    !     slip_particles_new, nparticle_slip, max_slip, dvec, gsvec, &
    !     lsvec, slip_likelihood_ls, slip_prior_ls, slip_weights, slip_mean, &
    !     slip_cov, slip_likelihood_ls_new, slip_prior_ls_new, &
    !     slip_st_rand, slip_particle_cur, slip_particle_cand, &
    !     slip_assigned_num, slip_id_start, slip_st_rand_ls, slip_metropolis_ls, &
    !     npath, nsar_index, nsar_total, sigma_sar_mat, gmat_arr, xmin, xmax, zmin, &
    !     fix_xbend, xbend, gsmc, qmean, qcov, qcov_inv, log_qcov_det)
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
        integer :: stride, offset_sample

        stride = 2
        offset_sample = 10

        slip_filepath = "output/more/mapslip_38400.dat"

        allocate (slip(2, nnode_total))
        allocate (particle_slip(ndim_slip))
        if (myid == 0) then
            print *, "fault_filepath: ", fault_filepath
            print *, "slip_filepath: ", slip_filepath
            allocate (fault_particles(ndim_fault, nparticle_fault))
            allocate (mean_slip(ndim_slip, nparticle_fault))
            allocate (tmp(ndim_fault + 1))
            open (10, file=fault_filepath, status='old')
            do k = 1, offset_sample
                read (10, *) tmp
            end do
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

            neglog = fault_calc_likelihood( &
                     particle, nplane, nxi_ls, neta_ls, nnode_total, ndof_total, ndof_index, ngnss, nobs, cny_fault, &
                     coor_fault, node_to_elem_val, node_to_elem_size, id_dof, luni, lmat, &
                     lmat_index, lmat_val, ltmat_index, ltmat_val, llmat, gmat, slip_dist, obs_points, &
                     obs_unitvec, obs_sigma, sigma2_full, alpha2_full, target_id_val, node_id_in_patch, &
                     xinode, etanode, uxinode, uetanode, r1vec, r2vec, nvec, response_dist, &
                     uobs, uret, slip_particles, slip_particles_new, &
                     nparticle_slip, max_slip, dvec, slip_likelihood_ls, slip_prior_ls, &
                     slip_weights, slip_mean, slip_cov, slip_likelihood_ls_new, &
                     slip_prior_ls_new, slip_assigned_num, slip_id_start, slip_st_rand_ls, &
                     slip_metropolis_ls, gsvec, lsvec, slip_particle_cur, &
                     slip_particle_cand, slip_st_rand, 0, "", &
                     npath, nsar_index, nsar_total, sigma_sar_mat, gmat_arr, xmin, xmax, zmin, fix_xbend, xbend)

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
                do inode = 1, nnode_total
                    do idirection = 1, 2
                        slip(idirection, inode) = 0d0
                    end do
                end do
                do idim = 1, ndim_slip
                    particle_slip(idim) = mean_slip(idim, iparticle)
                end do
                do idof = 1, ndof_total
                    inode = id_dof(idof)
                    slip(1, inode) = particle_slip(2*(idof - 1) + 1)
                    slip(2, inode) = particle_slip(2*(idof - 1) + 2)
                end do
                do inode = 1, nnode_total
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
        character(len=200) :: points_filepath, iplane_char
        double precision, allocatable :: fault_particles(:, :), mean_slip(:, :), &
            points(:, :)
        double precision :: r1, r2, xi, eta, uxi, ueta, x, y, z, &
            xf, yf, zf, strike, dip, lxi, leta, dip_rad, strike_rad, &
            theta(ndim_fault)

        integer :: ie, iparticle, idim, inode, ir1, ir2, jnode, idpoint, cnt, stride, k
        double precision :: pi = 4d0*atan(1d0)

        integer :: nr = 8, ir
        double precision, allocatable :: r1vec_cloud(:), r2vec_cloud(:)
        integer :: offset_patch, offset_sample

        allocate (r1vec_cloud(nr))
        allocate (r2vec_cloud(nr))
        do ir = 1, nr
            r1vec_cloud(ir) = -1d0 + 1d0/dble(nr) + (ir - 1)*2d0/dble(nr)
            r2vec_cloud(ir) = -1d0 + 1d0/dble(nr) + (ir - 1)*2d0/dble(nr)
        end do

        stride = 100
        offset_sample = 10
        points_filepath = "output/shuf/points_192000.dat"
        allocate (fault_particles(ndim_fault, nparticle_fault))
        ! x, y, z, uxi, ueta, unorm
        allocate (points(6, nparticle_fault*npatch_total*nr**2))
        allocate (mean_slip(2*nnode_total, nparticle_fault))
        allocate (tmp(ndim_fault + 1))
        allocate (particle(ndim_fault))
        allocate (svec(ndim_slip))
        open (10, file=fault_filepath, status='old')
        open (11, file=slip_filepath, status='old')
        do k = 1, offset_sample
            read (10, *) tmp
        end do
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
            slip_dist = reshape(mean_slip(:, iparticle), (/2, nnode_total/))
            do idim = 1, ndim_fault
                particle(idim) = fault_particles(idim, iparticle)
            end do

            call discretize_fault(particle, nplane, nxi_ls, neta_ls, cny_fault, coor_fault, &
                                  node_to_elem_val, node_to_elem_size, id_dof, xmin, xmax, zmin, fix_xbend, xbend)

            offset_patch = 0
            ! loop for patchs
            do iplane = 1, nplane
                nxi = nxi_ls(iplane)
                neta = neta_ls(iplane)
                xf = particle(8*(iplane - 1) + 1)
                yf = particle(8*(iplane - 1) + 2)
                zf = particle(8*(iplane - 1) + 3)
                strike = particle(8*(iplane - 1) + 4)
                strike_rad = strike*pi/180d0
                dip = particle(8*(iplane - 1) + 5)
                dip_rad = dip*pi/180d0
                lxi = particle(8*(iplane - 1) + 6)
                leta = particle(8*(iplane - 1) + 7)
                ! do ie = 1 + nxi*neta*(iplane - 1), nxi*neta*iplane
                do ie = offset_patch + 1, offset_patch + nxi*neta
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
                            idpoint = (cnt - 1)*npatch_total*nr*nr + (ie - 1)*nr*nr + (ir1 - 1)*nr + ir2
                            points(1, idpoint) = x
                            points(2, idpoint) = y
                            points(3, idpoint) = z
                            points(4, idpoint) = uxi
                            points(5, idpoint) = ueta
                            points(6, idpoint) = sqrt(uxi**2 + ueta**2)
                        end do
                    end do
                end do
                offset_patch = offset_patch + nxi*neta
            end do
            cnt = cnt + 1
        end do

        do iplane = 1, nplane
            write (iplane_char, "(I1)") iplane
            print *, trim(points_filepath)//trim(iplane_char)
            open (10 + iplane, file=trim(points_filepath)//trim(iplane_char), status='replace')
        end do
        do iparticle = 1, nparticle_fault
            offset_patch = 0
            do iplane = 1, nplane
                nxi = nxi_ls(iplane)
                neta = neta_ls(iplane)
                do i = (iparticle - 1)*npatch_total*nr**2 + offset_patch*nr**2 + 1, &
                    (iparticle - 1)*npatch_total*nr**2 + offset_patch*nr**2 + nxi*neta*nr**2
                    do j = 1, 6
                        write (10 + iplane, "(f12.5)", advance="no") points(j, i)
                    end do
                    write (10 + iplane, *)
                end do
                offset_patch = offset_patch + nxi*neta
            end do
            ! do i = , npatch_total*nr**2
            !     do j = 1, 6
            !         write (10, "(f12.5)", advance="no") points(j, i)
            !     end do
            !     write (10, *)
            ! end do
        end do
        do iplane = 1, nplane
            close (10 + iplane)
        end do
    end subroutine point_cloud

end program main
