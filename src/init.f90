module init
    implicit none
contains
    subroutine read_observation1(observation_path, &
                                 nobs)
        implicit none
        character(*), intent(in) :: observation_path
        integer, intent(inout) :: nobs
        character(len=200) :: buf
        integer :: ios
        open (10, file=observation_path)
        nobs = -1
        do
            read (10, *, iostat=ios) buf
            if (ios < 0) then
                exit
            end if
            nobs = nobs + 1
        end do
        close (10)
    end subroutine read_observation1

    subroutine read_observation2(observation_path, nobs, &
                                 obs_points, &
                                 obs_unitvec, obs_sigma, &
                                 dvec, nsar, ngnss)
        implicit none
        character(*), intent(in) :: observation_path
        integer, intent(in) :: nobs
        double precision, intent(inout) :: obs_points(:, :), &
            obs_unitvec(:, :), obs_sigma(:), dvec(:)
        integer, intent(inout) :: nsar, ngnss

        integer iobs
        double precision :: x, y, ex, ey, ez, dlos, sigma
        character(10) :: type
        character(200) :: buf
        open (10, file=observation_path)
        read (10, *) buf
        nsar = 0
        ngnss = 0
        ! for (iobs = 0 iobs < nobs iobs++) {
        do iobs = 1, nobs
            !     fscanf(fp, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%s\n", &x, &y, &ex, &ey, &ez,
            !            &dlos, &sigma, type)
            read (10, *) x, y, ex, ey, ez, dlos, sigma, type
            obs_points(1, iobs) = x
            obs_points(2, iobs) = y
            obs_unitvec(1, iobs) = ex
            obs_unitvec(2, iobs) = ey
            obs_unitvec(3, iobs) = ez
            dvec(iobs) = dlos
            obs_sigma(iobs) = sigma
            if (type == "sar") then
                nsar = nsar + 1
            end if
            if (type == "gnss") then
                ngnss = ngnss + 1
            end if
        end do
        ! GNSS observation have 3 direction components
        ngnss = ngnss/3
        close (10)
    end subroutine read_observation2

    subroutine discretize_fault(theta, nplane, nxi, neta, cny, coor, &
                                node_to_elem_val, node_to_elem_size, id_dof)
        implicit none
        double precision, intent(in) :: theta(:)
        integer, intent(in) :: nplane, nxi, neta
        double precision, intent(inout) :: coor(:, :)
        integer, intent(inout) :: cny(:, :), &
                                  node_to_elem_val(:, :), &
                                  node_to_elem_size(:), &
                                  id_dof(:)
        integer :: i, j, k, cnt, node_id, patch_id, iplane
        integer :: node1, node2, node3, node4
        double precision :: xi, eta, dxi, deta, lxi, leta

        cnt = 1
        do iplane = 1, nplane
            lxi = theta(iplane*8 - 2)
            leta = theta(iplane*8 - 1)
            ! length of a patch
            dxi = lxi/nxi
            deta = leta/neta

            ! coordinate of nodes
            do j = 1, neta + 1
                do i = 1, nxi + 1
                    node_id = i + (nxi + 1)*(j - 1) &
                              + (iplane - 1)*(nxi + 1)*(neta + 1)
                    xi = (i - 1)*dxi - lxi/2d0
                    eta = (j - 1)*deta - leta/2d0
                    coor(1, node_id) = xi
                    coor(2, node_id) = eta

                    node_to_elem_size(node_id) = 0
                    do k = 1, 4
                        node_to_elem_val(k, node_id) = -1
                    end do
                    ! no degree of freedom on the edge of the fault
                    if (i == 1 .or. i == nxi + 1 .or. &
                        j == 1 .or. j == neta + 1) then
                        cycle
                    end if
                    id_dof(cnt) = node_id
                    cnt = cnt + 1
                end do
            end do

            ! node id of patches
            do j = 1, neta
                do i = 1, nxi
                    patch_id = i + nxi*(j - 1) &
                               + (iplane - 1)*nxi*neta
                    node1 = i + (nxi + 1)*(j - 1) &
                            + (iplane - 1)*(nxi + 1)*(neta + 1)
                    node2 = node1 + 1
                    node3 = node2 + nxi + 1
                    node4 = node1 + nxi + 1
                    cny(1, patch_id) = node1
                    cny(2, patch_id) = node2
                    cny(3, patch_id) = node3
                    cny(4, patch_id) = node4
                    do k = 1, 4
                        node_id = cny(k, patch_id)
                        node_to_elem_size(node_id) = node_to_elem_size(node_id) + 1
                        node_to_elem_val(node_to_elem_size(node_id), node_id) &
                            = patch_id
                    end do
                end do
            end do
        end do
    end subroutine discretize_fault

    subroutine gen_laplacian(theta, nplane, nnode, nxi, neta, id_dof, ndof, luni, lmat)
        implicit none
        integer, intent(in) :: nplane, nnode, nxi, neta, id_dof(:), ndof
        double precision, intent(in) :: theta(:)
        double precision, intent(inout) ::  luni(:, :), lmat(:, :)
        double precision :: dxi, deta, lxi, leta
        integer :: iplane, inode, idof, jnode
        ! laplacian for single component
        do jnode = 1, nnode
            do inode = 1, nnode
                luni(inode, jnode) = 0d0
            end do
        end do

        do iplane = 1, nplane
            lxi = theta(iplane*8 - 2)
            leta = theta(iplane*8 - 1)
            dxi = lxi/nxi
            deta = leta/neta
            do inode = &
                1 + (nxi + 1)*(neta + 1)*(iplane - 1), (nxi + 1)*(neta + 1)*iplane
                if (mod(inode - 1, nxi + 1) == 0) then
                    luni(inode, inode + 1) = luni(inode, inode + 1) + 2d0/dxi**2d0
                    luni(inode, inode) = luni(inode, inode) - 2d0/dxi**2d0
                else if (mod(inode - 1, nxi + 1) == nxi) then
                    luni(inode, inode - 1) = luni(inode, inode - 1) + 2d0/dxi**2d0
                    luni(inode, inode) = luni(inode, inode) - 2d0/dxi**2d0
                else
                    luni(inode, inode - 1) = luni(inode, inode - 1) + 1d0/dxi**2d0
                    luni(inode, inode + 1) = luni(inode, inode + 1) + 1d0/dxi**2d0
                    luni(inode, inode) = luni(inode, inode) - 2d0/dxi**2d0
                end if

                if (mod((inode - 1)/(nxi + 1), neta + 1) == 0) then
                    luni(inode, inode + (nxi + 1)) = luni(inode, inode + (nxi + 1)) + 2d0/deta**2d0
                    luni(inode, inode) = luni(inode, inode) - 2d0/deta**2d0
                else if (mod((inode - 1)/(nxi + 1), neta + 1) == neta) then
                    luni(inode, inode - (nxi + 1)) = luni(inode, inode - (nxi + 1)) + 2d0/deta**2d0
                    luni(inode, inode) = luni(inode, inode) - 2d0/deta**2d0
                else
                    luni(inode, inode - (nxi + 1)) = luni(inode, inode - (nxi + 1)) + 1d0/deta**2d0
                    luni(inode, inode + (nxi + 1)) = luni(inode, inode + (nxi + 1)) + 1d0/deta**2d0
                    luni(inode, inode) = luni(inode, inode) - 2d0/deta**2d0
                end if
            end do
        end do

        ! laplacian for two components(u_xi, u_eta)
        do jnode = 1, 2*ndof
            do inode = 1, 2*nnode
                lmat(inode, jnode) = 0.
            end do
        end do

        do inode = 1, nnode
            do idof = 1, ndof
                jnode = id_dof(idof)
                lmat(2*(inode - 1) + 1, 2*(idof - 1) + 1) = luni(inode, jnode)
                lmat(2*(inode - 1) + 2, 2*(idof - 1) + 2) = luni(inode, jnode)
            end do
        end do
    end subroutine

    subroutine gen_sparse_lmat(lmat, lmat_index, lmat_val, ltmat_index, &
                               ltmat_val, nnode, ndof)
        implicit none
        integer, intent(in) :: nnode, ndof
        double precision, intent(in) :: lmat(:, :)
        integer, intent(inout) :: lmat_index(:, :), ltmat_index(:, :)
        double precision, intent(inout) :: lmat_val(:, :), ltmat_val(:, :)
        integer ::  i, j, cnt
        double precision :: val

        ! sparse matrix of L
        do i = 1, 2*nnode
            cnt = 0
            do j = 1, 2*ndof
                val = lmat(i, j)
                if (abs(val) > 1d-8) then
                    cnt = cnt + 1
                    lmat_index(cnt, i) = j
                    lmat_val(cnt, i) = val
                end if
            end do
            do while (cnt < 5)
                cnt = cnt + 1
                lmat_index(cnt, i) = 1
                lmat_val(cnt, i) = 0d0
            end do
        end do

        ! sparse matrix of L^T
        do i = 1, 2*ndof
            cnt = 0
            do j = 1, 2*nnode
                val = lmat(j, i)
                if (abs(val) > 1d-8) then
                    cnt = cnt + 1
                    ltmat_index(cnt, i) = j
                    ltmat_val(cnt, i) = val
                end if
            end do
            do while (cnt < 5)
                cnt = cnt + 1
                ltmat_index(cnt, i) = 1
                ltmat_val(cnt, i) = 0d0
            end do
        end do
    end subroutine gen_sparse_lmat

    subroutine calc_ll(llmat, lmat, nnode, ndof)
        implicit none
        integer, intent(in) :: nnode, ndof
        double precision, intent(in) :: lmat(:, :)
        double precision, intent(inout) :: llmat(:, :)
        double precision, allocatable :: tmp(:)
        integer :: n, m, i, j
        n = 2*nnode
        m = 2*ndof
        allocate (tmp(m))
        do j = 1, m
            do i = 1, m
                llmat(i, j) = 0d0
            end do
        end do
        call dgemm("t", "n", m, m, n, 1d0, lmat, n, &
                   lmat, n, 0d0, llmat, m)
    end subroutine calc_ll
end module init
