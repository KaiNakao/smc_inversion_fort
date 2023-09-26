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
    output_dir = "output_20000_20000_edge/"
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
    nparticle_slip = 20000
    nparticle_fault = 20000

    ! constrain max value for slip
    max_slip = 3d0

    ! dimension of fault parameter
    ndim_fault = 10
    ndim_slip = 2*ndof

    ! read observation data
    observation_path = "input/observation_toy2.csv"
    call read_observation1(observation_path, nobs)
    allocate (obs_points(2, nobs))
    allocate (obs_unitvec(3, nobs))
    allocate (obs_sigma(nobs))
    allocate (dvec(nobs))
    call read_observation2(observation_path, nobs, &
                           obs_points, &
                           obs_unitvec, obs_sigma, &
                           dvec, nsar, ngnss)

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
    ! particle = (/4.99413389000000, -9.81891359000000, -15.1300002400000, &
    !              0.872903704999999, 59.5438755200000, -6.545940000000001E-002, &
    !              0.350170955000000, -3.43126560500001, 31.6602526250000, &
    !              19.7599089800000/)
    ! print *, particle
    ! lxi = particle(9)
    ! leta = particle(10)
    ! call discretize_fault(lxi, leta, nxi, neta, cny_fault, coor_fault, &
    !                       node_to_elem_val, node_to_elem_size, id_dof)
    ! open (10, file="/home/nakao/smc_inversion/visualize/slip_mean.csv", &
    !       status='old')
    ! do i = 1, nnode
    !     read (10, *) slip(1, i), slip(2, i)
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
    ! open (10, file="/home/nakao/smc_inversion/visualize/dvec_est.dat", &
    !       status='replace')
    ! do i = 1, nobs
    !     write (10, *), gsvec(i)
    ! end do

    ! calculate likelihood for given fault
    allocate (particle(ndim_fault))
    ! allocate (tmp(ndim_fault + 1))
    ! do i = 1, ndim_fault
    !     particle(i) = 0d0
    ! end do
    ! do i = 1, ndim_fault + 1
    !     tmp(i) = 0d0
    ! end do
    ! open (10, file="output/14.csv", status='old')
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
    ! particle(:) = (/8.587707, -4.485232, -22.583390, 12.771108, 57.633782, &
    !                 -1.846341, 0.966436, 1.428768, 36.731036, 6.393948/)
    ! particle(:) = (/4.99413389000000, -9.81891359000000, -15.1300002400000, &
    !                 0.872903704999999, 59.5438755200000, -6.545940000000001E-002, &
    !                 0.350170955000000, -3.43126560500001, 31.6602526250000, &
    !                 19.7599089800000/)
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
    !          slip_particle_cand, slip_st_rand, 1, &
    !          "output/slip_from_mean_fault200000.dat")
    ! en_time = omp_get_wtime()
    ! print *, "etime: ", en_time - st_time
    ! print *, "neglog: ", neglog
    ! open (10, file="output/mean_faultsize.dat", status='replace')
    ! write (10, "(f12.5)") particle(9), particle(10)
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
    !          slip_particle_cand, slip_st_rand, 0, "")
    ! en_time = omp_get_wtime()
    ! print *, "etime: ", en_time - st_time
    ! print *, "neglog: ", neglog

    allocate (range(2, ndim_fault))
    range(:, :) = reshape((/-10., 10., -30., 0., -30., -1., -20., 20., 50., 90., &
                            -2., 2., -2., 2., -10., 2., 1., 50., 1., 50./), &
                          (/2, ndim_fault/))

    call fault_smc_exec( &
        output_dir, range, nparticle_fault, ndim_fault, &
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
    call mpi_finalize(ierr)

end program main
