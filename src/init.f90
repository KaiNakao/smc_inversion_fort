module init
    use type_mat
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
                                 dvec, nsar_ls, nsar_index, nsar_total, ngnss, npath)
        implicit none
        character(*), intent(in) :: observation_path
        integer, intent(in) :: nobs, npath
        double precision, intent(inout) :: obs_points(:, :), &
            obs_unitvec(:, :), obs_sigma(:), dvec(:)
        integer, intent(inout) :: nsar_ls(:), nsar_index(:), nsar_total, ngnss

        integer iobs, ipath
        double precision :: x, y, ex, ey, ez, dlos, sigma
        character(20) :: type, type_prev
        character(200) :: buf
        open (10, file=observation_path)
        read (10, *) buf

        nsar_total = 0
        ngnss = 0
        ipath = 1
        nsar_ls(ipath) = 0
        do iobs = 1, nobs
            read (10, *) x, y, ex, ey, ez, dlos, sigma, type
            obs_points(1, iobs) = x
            obs_points(2, iobs) = y
            obs_unitvec(1, iobs) = ex
            obs_unitvec(2, iobs) = ey
            obs_unitvec(3, iobs) = ez
            dvec(iobs) = dlos
            obs_sigma(iobs) = sigma
            if (iobs == 1) then
                type_prev = type
            end if
            if (type(1:5) == "InSAR") then
                if (type == type_prev) then
                    nsar_ls(ipath) = nsar_ls(ipath) + 1
                else
                    ipath = ipath + 1
                    nsar_ls(ipath) = 1
                end if
                type_prev = type
            end if
            if (type == "GNSS") then
                ngnss = ngnss + 1
            end if
        end do
        nsar_total = 0
        nsar_index(1) = 1
        do ipath = 1, npath
            nsar_total = nsar_total + nsar_ls(ipath)
            nsar_index(ipath + 1) = nsar_index(ipath) + nsar_ls(ipath)
        end do 
        ! GNSS observation have 3 direction components
        ngnss = ngnss/3
        close (10)
    end subroutine read_observation2

    subroutine calc_sigma_sar_mat(npath, nsar_ls, nsar_index, &
                                  obs_points, obs_sigma, sigma_sar_mat)
        implicit none
        integer, intent(in) :: npath, nsar_ls(:), nsar_index(:)
        double precision, intent(in) :: obs_points(:,:), obs_sigma(:)
        type(mat), intent(inout) :: sigma_sar_mat(:)
        integer :: ipath, iobs, jobs, nobs, iobs_start, lwork, info
        integer, allocatable :: ipiv(:)
        double precision, allocatable :: work(:)
        double precision :: xi, yi, xj, yj, dist

        do ipath = 1, npath
            nobs = nsar_ls(ipath)
            iobs_start = nsar_index(ipath)
            lwork = nobs
            allocate(ipiv(nobs))
            allocate(work(lwork))
            do iobs = 1, nobs
                xi = obs_points(1, iobs_start - 1 + iobs)
                yi = obs_points(2, iobs_start - 1 + iobs)
                do jobs = 1, nsar_ls(ipath)
                    xj = obs_points(1, iobs_start - 1 + jobs)
                    yj = obs_points(2, iobs_start - 1 + jobs)
                    dist = sqrt((xi - xj)**2 + (yi - yj)**2)
                    sigma_sar_mat(ipath)%body(iobs, jobs) = &
                        obs_sigma(nsar_index(ipath) + iobs) * &
                        obs_sigma(nsar_index(ipath) + jobs) * &
                        exp(-dist/1d1)
                end do
            end do
            call dgetrf(nobs, nobs, sigma_sar_mat(ipath)%body, &
                        nobs, ipiv, info)
            call dgetri(nobs, sigma_sar_mat(ipath)%body, &
                        nobs, ipiv, work, lwork, info)
            deallocate(ipiv)
            deallocate(work)
        end do
    end subroutine calc_sigma_sar_mat

    subroutine discretize_fault(theta, nplane, nxi_ls, neta_ls, cny, coor, &
                                node_to_elem_val, node_to_elem_size, id_dof)
        implicit none
        double precision, intent(in) :: theta(:)
        integer, intent(in) :: nplane, nxi_ls(:), neta_ls(:)
        double precision, intent(inout) :: coor(:, :)
        integer, intent(inout) :: cny(:, :), &
                                  node_to_elem_val(:, :), &
                                  node_to_elem_size(:), &
                                  id_dof(:)
        integer :: i, j, k, l, cnt, node_id, patch_id, iplane, inode
        integer :: offset_node, offset_patch
        integer :: nxi, neta
        integer :: node1, node2, node3, node4
        double precision :: xi, eta, dxi, deta, lxi, leta

        cnt = 1
        offset_node = 0
        offset_patch = 0
        do iplane = 1, nplane
            nxi = nxi_ls(iplane)
            neta = neta_ls(iplane)
            lxi = theta(iplane*8 - 2)
            leta = theta(iplane*8 - 1)
            ! length of a patch
            dxi = lxi/nxi
            deta = leta/neta

            ! coordinate of nodes
            do j = 1, neta + 1
                do i = 1, nxi + 1
                    ! node_id = i + (nxi + 1)*(j - 1) &
                    !           + (iplane - 1)*(nxi + 1)*(neta + 1) + offset
                    node_id = i + (nxi + 1)*(j - 1) + offset_node
                    xi = (i - 1)*dxi - lxi/2d0
                    eta = (j - 1)*deta - leta/2d0
                    coor(1, node_id) = xi
                    coor(2, node_id) = eta

                    node_to_elem_size(node_id) = 0
                    do k = 1, 4
                        node_to_elem_val(k, node_id) = -1
                    end do
                    ! no degree of freedom on the edge of the fault
                    ! if (iplane == 1) then
                    !     if (i == 1 .or. j == 1) then
                    !         cycle
                    !     end if
                    ! else if (iplane == nplane) then
                    !     if (i == nxi + 1 .or. j == 1) then
                    !         cycle
                    !     end if
                    ! else
                    !     if (j == 1) then
                    !         cycle
                    !     end if
                    ! end if
                    if (j == 1) then
                        cycle
                    end if
                    if (iplane == 1) then
                        if (i == 1) then
                            cycle
                        end if
                    end if
                    if (iplane == nplane) then
                        if (i == nxi + 1) then
                            cycle
                        end if
                    end if
                    id_dof(cnt) = node_id
                    cnt = cnt + 1
                end do
            end do

            ! node id of patches
            do j = 1, neta
                do i = 1, nxi
                    ! patch_id = i + nxi*(j - 1) &
                    !            + (iplane - 1)*nxi*neta
                    patch_id = offset_patch + i + nxi*(j - 1)
                    ! node1 = i + (nxi + 1)*(j - 1) &
                    !         + (iplane - 1)*(nxi + 1)*(neta + 1)
                    node1 = i + (nxi + 1)*(j - 1) &
                            + offset_node
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
            offset_patch = offset_patch + nxi * neta
            offset_node = offset_node + (nxi + 1) * (neta + 1)
        end do

        ! offset_node = 0
        ! do iplane = 1, nplane
        !     nxi = nxi_ls(iplane)
        !     neta = neta_ls(iplane)
        !     do inode = offset_node + 1, offset_node + (nxi + 1) * (neta + 1)
        !         k = mod((inode - 1 - offset_node), (nxi + 1)) + 1
        !         l = (inode - 1 - offset_node) / (nxi + 1) + 1
        !         print *, coor(1, inode), coor(2, inode), k, l
        !     end do
        !     offset_node = offset_node + (nxi + 1) * (neta + 1)
        ! end do
        ! stop
    end subroutine discretize_fault

    subroutine gen_laplacian(theta, nplane, nnode_total, nxi_ls, neta_ls, id_dof, ndof_total, luni, lmat)
        implicit none
        integer, intent(in) :: nplane, nnode_total, nxi_ls(:), neta_ls(:), id_dof(:), ndof_total
        double precision, intent(in) :: theta(:)
        double precision, intent(inout) ::  luni(:, :), lmat(:, :)
        double precision :: dxi, deta, lxi, leta, dcross
        integer :: iplane, inode, idof, jnode, nxi, neta, offset, k, l
        ! laplacian for single component
        do jnode = 1, nnode_total
            do inode = 1, nnode_total
                luni(inode, jnode) = 0d0
            end do
        end do

        offset = 0
        do iplane = 1, nplane
            nxi = nxi_ls(iplane)
            neta = neta_ls(iplane)
            lxi = theta(iplane*8 - 2)
            leta = theta(iplane*8 - 1)
            dxi = lxi/nxi
            deta = leta/neta
            dcross = sqrt(dxi**2 + deta**2)
            ! do inode = &
            !     1 + (nxi + 1)*(neta + 1)*(iplane - 1), (nxi + 1)*(neta + 1)*iplane
            do inode = offset + 1, offset + (nxi + 1) * (neta + 1)
                k = mod((inode - 1 - offset), (nxi + 1)) + 1
                l = (inode - 1 - offset) / (nxi + 1) + 1

                luni(inode, inode) = luni(inode, inode) - 2d0/dxi**2
                if (k == nxi + 1) then
                    luni(inode, inode - 1) = luni(inode, inode - 1) + 1d0/dxi**2
                else 
                    luni(inode, inode + 1) = luni(inode, inode + 1) + 1d0/dxi**2
                end if

                if (k == 1) then
                    luni(inode, inode + 1) = luni(inode, inode + 1) + 1d0/dxi**2
                else 
                    luni(inode, inode - 1) = luni(inode, inode - 1) + 1d0/dxi**2
                end if

                if (l == neta + 1) then
                    luni(inode, inode - (nxi + 1)) = luni(inode, inode - (nxi + 1)) + 1d0/deta**2
                else 
                    luni(inode, inode + (nxi + 1)) = luni(inode, inode + (nxi + 1)) + 1d0/deta**2
                end if

                if (l == 1) then
                    luni(inode, inode + (nxi + 1)) = luni(inode, inode + (nxi + 1)) + 1d0/deta**2
                else 
                    luni(inode, inode - (nxi + 1)) = luni(inode, inode - (nxi + 1)) + 1d0/deta**2
                end if

                if (k == nxi + 1 .and. l == neta + 1) then
                    luni(inode, inode - (nxi + 1) - 1) = luni(inode, inode - (nxi + 1) - 1) + 1d0/(4d0 * dxi * deta)
                else if (k == nxi + 1) then
                    luni(inode, inode + (nxi + 1) - 1) = luni(inode, inode + (nxi + 1) - 1) + 1d0/(4d0 * dxi * deta)
                else if (l == neta + 1) then
                    luni(inode, inode - (nxi + 1) + 1) = luni(inode, inode - (nxi + 1) + 1) + 1d0/(4d0 * dxi * deta)
                else
                    luni(inode, inode + (nxi + 1) + 1) = luni(inode, inode + (nxi + 1) + 1) + 1d0/(4d0 * dxi * deta)
                end if

                if (k == 1 .and. l == neta + 1) then
                    luni(inode, inode - (nxi + 1) + 1) = luni(inode, inode - (nxi + 1) + 1) - 1d0/(4d0 * dxi * deta)
                else if (k == 1) then
                    luni(inode, inode + (nxi + 1) + 1) = luni(inode, inode + (nxi + 1) + 1) - 1d0/(4d0 * dxi * deta)
                else if (l == neta + 1) then
                    luni(inode, inode - (nxi + 1) - 1) = luni(inode, inode - (nxi + 1) - 1) - 1d0/(4d0 * dxi * deta)
                else
                    luni(inode, inode + (nxi + 1) - 1) = luni(inode, inode + (nxi + 1) - 1) - 1d0/(4d0 * dxi * deta)
                end if

                if (k == 1 .and. l == 1) then
                    luni(inode, inode + (nxi + 1) + 1) = luni(inode, inode + (nxi + 1) + 1) + 1d0/(4d0 * dxi * deta)
                else if (k == 1) then
                    luni(inode, inode - (nxi + 1) + 1) = luni(inode, inode - (nxi + 1) + 1) + 1d0/(4d0 * dxi * deta)
                else if (l == 1) then
                    luni(inode, inode + (nxi + 1) - 1) = luni(inode, inode + (nxi + 1) - 1) + 1d0/(4d0 * dxi * deta)
                else
                    luni(inode, inode - (nxi + 1) - 1) = luni(inode, inode - (nxi + 1) - 1) + 1d0/(4d0 * dxi * deta)
                end if

                if (k == nxi + 1 .and. l == 1) then
                    luni(inode, inode + (nxi + 1) - 1) = luni(inode, inode + (nxi + 1) - 1) - 1d0/(4d0 * dxi * deta)
                else if (k == nxi + 1) then
                    luni(inode, inode - (nxi + 1) - 1) = luni(inode, inode - (nxi + 1) - 1) - 1d0/(4d0 * dxi * deta)
                else if (l == 1) then
                    luni(inode, inode + (nxi + 1) + 1) = luni(inode, inode + (nxi + 1) + 1) - 1d0/(4d0 * dxi * deta)
                else
                    luni(inode, inode - (nxi + 1) + 1) = luni(inode, inode - (nxi + 1) + 1) - 1d0/(4d0 * dxi * deta)
                end if
                ! if (k == 1) then
                !     luni(inode, inode + 1) = luni(inode, inode + 1) + 2d0/dxi**2
                ! else if (k == nxi + 1) then
                !     luni(inode, inode - 1) = luni(inode, inode - 1) + 2d0/dxi**2
                ! else
                !     luni(inode, inode + 1) = luni(inode, inode + 1) + 1d0/dxi**2
                !     luni(inode, inode - 1) = luni(inode, inode - 1) + 1d0/dxi**2
                ! end if

                ! luni(inode, inode) = luni(inode, inode) - 2d0/deta**2
                ! if (l == 1) then
                !     luni(inode, inode + (nxi + 1)) = luni(inode, inode + (nxi + 1)) + 2d0/deta**2
                ! else if (l == neta + 1) then
                !     luni(inode, inode - (nxi + 1)) = luni(inode, inode - (nxi + 1)) + 2d0/deta**2
                ! else
                !     luni(inode, inode + (nxi + 1)) = luni(inode, inode + (nxi + 1)) + 1d0/deta**2
                !     luni(inode, inode - (nxi + 1)) = luni(inode, inode - (nxi + 1)) + 1d0/deta**2
                ! end if

                ! luni(inode, inode) = luni(inode, inode) - 2d0/dcross**2
                ! if ((k == 1 .and. l < neta + 1) .or. (k < nxi + 1 .and. l == 1)) then
                !     luni(inode, inode + (nxi + 1) + 1) = luni(inode, inode + (nxi + 1) + 1) + 2d0/dcross**2
                ! else if ((k == nxi + 1 .and. l > 1) .or. (k > 1 .and. l == neta + 1)) then
                !     luni(inode, inode - (nxi + 1) + 1) = luni(inode, inode - (nxi + 1) + 1) + 2d0/dcross**2
                ! else if (k > 1 .and. k < nxi + 1 .and. l > 1 .and. l < neta + 1) then
                !     luni(inode, inode + (nxi + 1) + 1) = luni(inode, inode + (nxi + 1) + 1) + 1d0/dcross**2
                !     luni(inode, inode - (nxi + 1) + 1) = luni(inode, inode - (nxi + 1) + 1) + 1d0/dcross**2
                ! end if

                ! luni(inode, inode) = luni(inode, inode) - 2d0/dcross**2
                ! if ((k == 1 .and. l > 1) .or. (k < nxi + 1 .and. l == neta + 1)) then
                !     luni(inode, inode - (nxi + 1) + 1) = luni(inode, inode - (nxi + 1) + 1) + 2d0/dcross**2
                ! else if ((k == nxi + 1 .and. l < neta + 1) .or. (k > 1 .and. l == 1)) then
                !     luni(inode, inode + (nxi + 1) - 1) = luni(inode, inode + (nxi + 1) - 1) + 2d0/dcross**2
                ! else if (k > 1 .and. k < nxi + 1 .and. l > 1 .and. l < neta + 1) then
                !     luni(inode, inode - (nxi + 1) + 1) = luni(inode, inode - (nxi + 1) + 1) + 1d0/dcross**2
                !     luni(inode, inode + (nxi + 1) - 1) = luni(inode, inode + (nxi + 1) - 1) + 1d0/dcross**2
                ! end if
                ! if (mod(inode - 1 - offset, nxi + 1) == 0) then
                !     luni(inode, inode + 1) = luni(inode, inode + 1) + 2d0/dxi**2d0
                !     luni(inode, inode) = luni(inode, inode) - 2d0/dxi**2d0
                ! else if (mod(inode - 1 - offset, nxi + 1) == nxi) then
                !     luni(inode, inode - 1) = luni(inode, inode - 1) + 2d0/dxi**2d0
                !     luni(inode, inode) = luni(inode, inode) - 2d0/dxi**2d0
                ! else
                !     luni(inode, inode - 1) = luni(inode, inode - 1) + 1d0/dxi**2d0
                !     luni(inode, inode + 1) = luni(inode, inode + 1) + 1d0/dxi**2d0
                !     luni(inode, inode) = luni(inode, inode) - 2d0/dxi**2d0
                ! end if

                ! if (mod((inode - 1 - offset)/(nxi + 1), neta + 1) == 0) then
                !     luni(inode, inode + (nxi + 1)) = luni(inode, inode + (nxi + 1)) + 2d0/deta**2d0
                !     luni(inode, inode) = luni(inode, inode) - 2d0/deta**2d0
                ! else if (mod((inode - 1 - offset)/(nxi + 1), neta + 1) == neta) then
                !     luni(inode, inode - (nxi + 1)) = luni(inode, inode - (nxi + 1)) + 2d0/deta**2d0
                !     luni(inode, inode) = luni(inode, inode) - 2d0/deta**2d0
                ! else
                !     luni(inode, inode - (nxi + 1)) = luni(inode, inode - (nxi + 1)) + 1d0/deta**2d0
                !     luni(inode, inode + (nxi + 1)) = luni(inode, inode + (nxi + 1)) + 1d0/deta**2d0
                !     luni(inode, inode) = luni(inode, inode) - 2d0/deta**2d0
                ! end if
            end do
            offset = offset + (nxi + 1) * (neta + 1)
        end do

        ! laplacian for two components(u_xi, u_eta)
        do jnode = 1, 2*ndof_total
            do inode = 1, 2*nnode_total
                lmat(inode, jnode) = 0d0
            end do
        end do

        do inode = 1, nnode_total
            do idof = 1, ndof_total
                jnode = id_dof(idof)
                lmat(2*(inode - 1) + 1, 2*(idof - 1) + 1) = luni(inode, jnode)
                lmat(2*(inode - 1) + 2, 2*(idof - 1) + 2) = luni(inode, jnode)
            end do
        end do
    end subroutine

    subroutine gen_sparse_lmat(lmat, lmat_index, lmat_val, ltmat_index, &
                               ltmat_val, nnode_total, ndof_total)
        implicit none
        integer, intent(in) :: nnode_total, ndof_total
        double precision, intent(in) :: lmat(:, :)
        integer, intent(inout) :: lmat_index(:, :), ltmat_index(:, :)
        double precision, intent(inout) :: lmat_val(:, :), ltmat_val(:, :)
        integer ::  i, j, cnt
        double precision :: val

        ! sparse matrix of L
        do i = 1, 2*nnode_total
            cnt = 0
            do j = 1, 2*ndof_total
                val = lmat(i, j)
                if (abs(val) > 1d-8) then
                    cnt = cnt + 1
                    lmat_index(cnt, i) = j
                    lmat_val(cnt, i) = val
                end if
            end do
            do while (cnt < 9)
                cnt = cnt + 1
                lmat_index(cnt, i) = 1
                lmat_val(cnt, i) = 0d0
            end do
        end do

        ! sparse matrix of L^T
        do i = 1, 2*ndof_total
            cnt = 0
            do j = 1, 2*nnode_total
                val = lmat(j, i)
                if (abs(val) > 1d-8) then
                    cnt = cnt + 1
                    ltmat_index(cnt, i) = j
                    ltmat_val(cnt, i) = val
                end if
            end do
            do while (cnt < 9)
                cnt = cnt + 1
                ltmat_index(cnt, i) = 1
                ltmat_val(cnt, i) = 0d0
            end do
        end do
    end subroutine gen_sparse_lmat

    subroutine calc_ll(llmat, lmat, nnode_total, ndof_total)
        implicit none
        integer, intent(in) :: nnode_total, ndof_total
        double precision, intent(in) :: lmat(:, :)
        double precision, intent(inout) :: llmat(:, :)
        double precision, allocatable :: tmp(:)
        integer :: n, m, i, j
        n = 2*nnode_total
        m = 2*ndof_total
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
