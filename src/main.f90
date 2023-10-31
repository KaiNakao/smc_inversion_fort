program main
    use mpi
    use omp_lib
    use init
    use gfunc
    use smc_fault
    implicit none
    integer :: i, j
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
    double precision :: xf, yf, zf, strike, dip, log_sigma_sar2, &
        log_sigma_gnss2, log_alpha2, lxi, leta

    ! laplacian
    double precision, allocatable :: luni(:, :), lmat(:, :), &
        lmat_val(:, :), llmat(:, :)
    integer, allocatable :: lmat_index(:, :)

    ! SMC for slip
    double precision, allocatable :: slip_particles(:, :), &
        slip_particles_new(:, :), slip_likelihood_ls(:), slip_prior_ls(:), &
        slip_weights(:), slip_mean(:), slip_cov(:, :), slip_likelihood_ls_new(:), &
        slip_prior_ls_new(:), slip_st_rand(:), slip_st_rand_ls(:, :), &
        slip_metropolis_ls(:), sigma2_full(:), gsvec(:), lsvec(:), &
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

    ! set MPI
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world, myid, ierr)
    call mpi_comm_size(mpi_comm_world, numprocs, ierr)

    ! random seed
    call random_seed(size=seedsize)
    allocate (seed(seedsize))
    do i = 1, seedsize
        seed(i) = 2
    end do
    call random_seed(put=seed)

    ! generate output directory
    output_dir = "output_200000_192000/"

    command = "mkdir -p "
    command = trim(command)//" "//output_dir
    call system(command)

    ! number of patches
    nxi = 6
    neta = 6

    ! number of patches, nodes, dof
    npatch = nxi*neta
    nnode = (nxi + 1)*(neta + 1)
    ndof = (nxi - 1)*(neta - 1)

    ! number of samples
    nparticle_slip = 200000
    nparticle_fault = 20000

    ! constrain max value for slip
    max_slip = 3d0

    ! dimension of fault parameter
    ndim_fault = 10
    ndim_slip = 2*ndof

    ! read observation data
    ! synthetic test
    ! observation_path = "input/observation_toy2.csv"
    ! real observation data
    observation_path = "input/observation_with_gnss_reduced.csv"
    if (myid == 0) then
        call read_observation1(observation_path, nobs)
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
    allocate (gmat(nobs, 2*ndof))
    allocate (slip_dist(2, nnode))
    allocate (response_dist(3, nobs))

    ! laplacian
    allocate (luni(nnode, nnode))
    allocate (lmat(2*nnode, 2*ndof))
    allocate (lmat_index(5, 2*nnode))
    allocate (lmat_val(5, 2*nnode))
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
    allocate (gsvec(nobs))
    allocate (lsvec(2*nnode))
    allocate (slip_particle_cur(ndim_slip))
    allocate (slip_particle_cand(ndim_slip))

    ! ! calculate diplacement for fixed fault and slip distribution
    ! allocate (slip(2, nnode))
    ! allocate (svec(2*ndof))

    ! allocate (particle(ndim_fault))
    ! open (10, file="visualize/mean_fault.dat", status='old')
    ! do i = 1, ndim_fault
    !     read (10, *) particle(i)
    !     print *, i, particle(i)
    ! end do
    ! close (10)

    ! print *, particle
    ! lxi = particle(9)
    ! leta = particle(10)
    ! call discretize_fault(lxi, leta, nxi, neta, cny_fault, coor_fault, &
    !                       node_to_elem_val, node_to_elem_size, id_dof)
    ! open (10, file="visualize/mean_slip.dat", status='old')
    ! do i = 1, nnode
    !     read (10, *) slip(1, i), slip(2, i)
    !     print *, i, slip(1, i), slip(2, i)
    ! end do
    ! do i = 1, ndof
    !     svec(2*i - 1) = slip(1, id_dof(i))
    !     svec(2*i) = slip(2, id_dof(i))
    ! end do
    ! close (10)

    ! print *, (svec)
    ! xf = particle(1)
    ! yf = particle(2)
    ! zf = particle(3)
    ! strike = particle(4)
    ! dip = particle(5)
    ! log_sigma_sar2 = particle(6)
    ! log_sigma_gnss2 = particle(7)
    ! log_alpha2 = particle(8)
    ! call calc_greens_func(gmat, slip_dist, cny_fault, coor_fault, &
    !                       obs_points, obs_unitvec, leta, xf, yf, zf, strike, dip, &
    !                       node_to_elem_val, node_to_elem_size, id_dof, nsar, ngnss, &
    !                       nobs, nnode, ndof, target_id_val, node_id_in_patch, xinode, &
    !                       etanode, uxinode, uetanode, r1vec, r2vec, nvec, response_dist, &
    !                       uobs, uret)

    ! call dgemv('n', nobs, 2*ndof, 1d0, gmat, &
    !            nobs, svec, 1, 0d0, gsvec, 1)
    ! open (10, file="visualize/dvec_est.dat", status='replace')
    ! do i = 1, nobs
    !     write (10, *), gsvec(i)
    ! end do

    ! ! calculate likelihood for given fault
    ! allocate (particle(ndim_fault))
    ! allocate (tmp(ndim_fault + 1))
    ! do i = 1, ndim_fault
    !     particle(i) = 0d0
    ! end do
    ! do i = 1, ndim_fault + 1
    !     tmp(i) = 0d0
    ! end do
    ! open (10, file="output_200000_192000_iburi/23.csv", status='old')
    ! do i = 1, nparticle_fault
    !     read (10, *) tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), &
    !         tmp(6), tmp(7), tmp(8), tmp(9), tmp(10), tmp(11)
    !     do j = 1, ndim_fault
    !         particle(j) = particle(j) + tmp(j)
    !     end do
    ! end do
    ! close (10)
    ! do j = 1, ndim_fault
    !     particle(j) = particle(j)/nparticle_fault
    ! end do
    ! print *, particle
    ! open (10, file="visualize/mean_fault.dat", status="replace")
    ! do i = 1, ndim_fault
    !     write (10, *) particle(i)
    ! end do
    ! close (10)

    ! st_time = omp_get_wtime()
    ! neglog = fault_calc_likelihood( &
    !          particle, nxi, neta, nnode, ndof, nsar, ngnss, nobs, cny_fault, &
    !          coor_fault, node_to_elem_val, node_to_elem_size, id_dof, luni, lmat, &
    !          lmat_index, lmat_val, llmat, gmat, slip_dist, obs_points, &
    !          obs_unitvec, obs_sigma, sigma2_full, target_id_val, node_id_in_patch, &
    !          xinode, etanode, uxinode, uetanode, r1vec, r2vec, nvec, response_dist, &
    !          uobs, uret, slip_particles, slip_particles_new, &
    !          nparticle_slip, max_slip, dvec, slip_likelihood_ls, slip_prior_ls, &
    !          slip_weights, slip_mean, slip_cov, slip_likelihood_ls_new, &
    !          slip_prior_ls_new, slip_assigned_num, slip_id_start, slip_st_rand_ls, &
    !          slip_metropolis_ls, gsvec, lsvec, slip_particle_cur, &
    !          slip_particle_cand, slip_st_rand, 1, "visualize/slip_from_mean_fault.dat")
    ! en_time = omp_get_wtime()
    ! print *, "etime: ", en_time - st_time
    ! print *, "neglog: ", neglog

    call calc_slip_map("var_faultsize_max3/23_converted.csv")
    ! call point_cloud("var_faultsize_max3/23_converted.csv", &
    !                  "var_faultsize_max3/mapslip.dat")

    ! ! smc for fault
    ! allocate (range(2, ndim_fault))
    ! ! range of uniform prior distribution P(theta)
    ! ! synthetic test
    ! ! range(:, :) = reshape((/-5., 15., -15., 15., -39., -10., -20., 20., 50., 90., &
    ! !                         -2., 2., -2., 2., -10., 2., 1., 50., 1., 50./), &
    ! !                       (/2, ndim_fault/))
    ! ! real observation data
    ! range(:, :) = reshape((/-10, 10, -30, 0, -30, -1, -20, 20, 50, 90, &
    !                         -2, 2, -2, 2, -10, 2, 1, 50, 1, 50/), &
    !                       (/2, ndim_fault/))
    ! call fault_smc_exec( &
    !     output_dir, range, nparticle_fault, ndim_fault, &
    !     myid, numprocs, nxi, neta, nnode, ndof, nsar, ngnss, nobs, cny_fault, &
    !     coor_fault, node_to_elem_val, node_to_elem_size, id_dof, luni, lmat, &
    !     lmat_index, lmat_val, llmat, gmat, slip_dist, obs_points, &
    !     obs_unitvec, obs_sigma, sigma2_full, target_id_val, node_id_in_patch, &
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

        slip_filepath = "var_faultsize_max3/mapslip.dat"
        allocate (slip(2, nnode))
        allocate (particle_slip(ndim_slip))
        if (myid == 0) then
            allocate (fault_particles(ndim_fault, nparticle_fault))
            allocate (mean_slip(ndim_slip, nparticle_fault))
            allocate (tmp(ndim_fault + 1))
            open (10, file=fault_filepath, status='old')
            do i = 1, nparticle_fault
                read (10, *) tmp
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
                     particle, nxi, neta, nnode, ndof, nsar, ngnss, nobs, cny_fault, &
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
        double precision :: M_PI, dip_rad, strike_rad, r1, r2, xi, eta, uxi, ueta, &
            x, y, z
        integer :: ie, iparticle, idim, inode, ir1, ir2, jnode, idpoint, cnt, stride

        stride = 5
        points_filepath = "var_faultsize_max3/points.dat"
        allocate (fault_particles(ndim_fault, nparticle_fault))
        ! x, y, z, uxi, ueta, unorm
        allocate (points(6, nparticle_fault*npatch*4))
        allocate (mean_slip(2*nnode, nparticle_fault))
        allocate (tmp(ndim_fault + 1))
        allocate (particle(ndim_fault))
        allocate (svec(ndim_slip))
        open (10, file=fault_filepath, status='old')
        open (11, file=slip_filepath, status='old')
        do i = 1, nparticle_fault
            read (10, *) tmp
            read (11, *) mean_slip(:, i)
            do j = 1, ndim_fault
                fault_particles(j, i) = tmp(j)
            end do
        end do
        close (10)
        close (11)

        cnt = 1
        do iparticle = 1, nparticle_fault, stride
            slip_dist = reshape(mean_slip(:, iparticle), (/2, nnode/))
            do idim = 1, ndim_fault
                particle(idim) = fault_particles(idim, iparticle)
            end do
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

            call discretize_fault(lxi, leta, nxi, neta, cny_fault, coor_fault, &
                                  node_to_elem_val, node_to_elem_size, id_dof)

            M_PI = 2d0*asin(1d0)
            dip_rad = dip/180d0*M_PI
            strike_rad = strike/180d0*M_PI
            r1vec(1) = -0.5d0
            r1vec(2) = 0.5d0
            r2vec(1) = -0.5d0
            r2vec(2) = 0.5d0
            ! loop for patchs
            do ie = 1, npatch
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
                do ir1 = 1, 2
                    r1 = r1vec(ir1)
                    do ir2 = 1, 2
                        r2 = r2vec(ir2)
                        !   interpolate(xi, eta, uxi, ueta) by shape functions
                        nvec(1) = (1d0 - r1)*(1d0 - r2)/4d0
                        nvec(2) = (1d0 + r1)*(1d0 - r2)/4d0
                        nvec(3) = (1d0 + r1)*(1d0 + r2)/4d0
                        nvec(4) = (1d0 - r1)*(1d0 + r2)/4d0
                        xi = 0d0
                        eta = 0d0
                        uxi = 0d0
                        ueta = 0d0
                        !   for(inode=0 inode < 4 inode + +) {
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
            cnt = cnt + 1
        end do

        open (10, file=points_filepath, status='replace')
        do i = 1, nparticle_fault*npatch*4/stride
            do j = 1, 6
                write (10, "(f12.5)", advance="no") points(j, i)
            end do
            write (10, *)
        end do
        close (10)
    end subroutine point_cloud

    ! subroutine point_cloud(fault_filepath, slip_filepath)
    !     implicit none
    !     character(*), intent(in) :: fault_filepath, slip_filepath
    !     character(len=200) :: points_filepath
    !     double precision, allocatable :: fault_particles(:, :), mean_slip(:, :), &
    !         points(:, :)
    !     double precision :: M_PI, dip_rad, strike_rad, r1, r2, xi, eta, uxi, ueta, &
    !         x, y, z
    !     integer :: ie, iparticle, idim, inode, ir1, ir2, jnode, idpoint, cnt, stride

    !     stride = 10
    !     points_filepath = "output_200000_192000_iburi/points.dat"
    !     allocate (fault_particles(ndim_fault, nparticle_fault))
    !     ! x, y, z, uxi, ueta, unorm
    !     allocate (points(6, nparticle_fault*nnode))
    !     allocate (mean_slip(2*nnode, nparticle_fault))
    !     allocate (tmp(ndim_fault + 1))
    !     allocate (particle(ndim_fault))
    !     open (10, file=fault_filepath, status='old')
    !     open (11, file=slip_filepath, status='old')
    !     do i = 1, nparticle_fault
    !         read (10, *) tmp
    !         do j = 1, ndim_fault
    !             fault_particles(j, i) = tmp(j)
    !         end do
    !         read (11, *) mean_slip(:, i)
    !     end do
    !     close (10)
    !     close (11)

    !     cnt = 1
    !     do iparticle = 1, nparticle_fault, stride
    !         slip_dist = reshape(mean_slip(:, iparticle), (/2, nnode/))
    !         do idim = 1, ndim_fault
    !             particle(idim) = fault_particles(idim, iparticle)
    !         end do
    !         xf = particle(1)
    !         yf = particle(2)
    !         zf = particle(3)
    !         strike = particle(4)
    !         dip = particle(5)
    !         log_sigma_sar2 = particle(6)
    !         log_sigma_gnss2 = particle(7)
    !         log_alpha2 = particle(8)
    !         lxi = particle(9)
    !         leta = particle(10)

    !         call discretize_fault(lxi, leta, nxi, neta, cny_fault, coor_fault, &
    !                               node_to_elem_val, node_to_elem_size, id_dof)

    !         M_PI = 2d0*asin(1d0)
    !         dip_rad = dip/180d0*M_PI
    !         strike_rad = strike/180d0*M_PI
    !         ! loop for nodes
    !         do inode = 1, nnode
    !             xi = coor_fault(1, inode)
    !             eta = coor_fault(2, inode)
    !             uxi = slip_dist(1, inode)
    !             ueta = slip_dist(2, inode)
    !             x = xf - eta*cos(dip_rad)*cos(strike_rad) + &
    !                 xi*sin(strike_rad)
    !             y = yf + eta*cos(dip_rad)*sin(strike_rad) + &
    !                 xi*cos(strike_rad)
    !             z = zf + eta*sin(dip_rad)
    !             idpoint = (cnt - 1)*nnode + inode
    !             points(1, idpoint) = x
    !             points(2, idpoint) = y
    !             points(3, idpoint) = z
    !             points(4, idpoint) = uxi
    !             points(5, idpoint) = ueta
    !             points(6, idpoint) = sqrt(uxi**2 + ueta**2)
    !         end do
    !         cnt = cnt + 1
    !     end do

    !     open (10, file=points_filepath, status='replace')
    !     do i = 1, nparticle_fault*nnode/stride
    !         do j = 1, 6
    !             write (10, "(f12.5)", advance="no") points(j, i)
    !         end do
    !         write (10, *)
    !     end do
    !     close (10)
    ! end subroutine point_cloud
end program main
