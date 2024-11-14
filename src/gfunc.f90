module gfunc
    use init
    use type_mat
    implicit none
contains
    subroutine gen_unit_slip(inode, idirection, slip_dist)
        implicit none
        integer, intent(in) :: inode, idirection
        double precision, intent(inout) :: slip_dist(:, :)
        ! slip value at nodes(u_xi, u_eta components)
        ! in Iburi case, "the value of u_xi"is expected to be negative.
        ! but it need to be handled with positivity constraints
        if (idirection == 1) then
            slip_dist(idirection, inode) = -1d0
        end if
        if (idirection == 2) then
            slip_dist(idirection, inode) = 1d0
        end if
    end subroutine gen_unit_slip

    subroutine del_unit_slip(inode, idirection, slip_dist)
        implicit none
        integer, intent(in) :: inode, idirection
        double precision, intent(inout) :: slip_dist(:, :)
        slip_dist(idirection, inode) = 0d0
    end subroutine del_unit_slip

    ! subroutine call_dc3d0(xsource, ysource, zsource, xobs, yobs, uxi, ueta, &
    !                       dip, area, strike, uret)
    !     implicit none
    !     double precision, intent(in) :: xsource, ysource, zsource, xobs, yobs, &
    !         uxi, ueta, area, strike
    !     double precision, intent(inout) :: uret(:)
    !     real, intent(inout) :: dip
    !     double precision :: M_PI, strike_rad, ux_rot, uy_rot
    !     real :: x, y, alpha, z, depth, pot1, pot2, pot3, pot4, &
    !             ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz
    !     integer :: iret
    !     M_PI = 2d0*asin(1d0)
    !     strike_rad = strike/180d0*M_PI
    !     ! in DC3D, x axis is the strike direction
    !     ! rotation
    !     x = real(sin(strike_rad)*(xobs - xsource) + cos(strike_rad)*(yobs - ysource))
    !     y = real(-cos(strike_rad)*(xobs - xsource) + sin(strike_rad)*(yobs - ysource))
    !     alpha = real(2d0/3d0)
    !     z = 0d0
    !     depth = -real(zsource)
    !     pot1 = real(uxi*area)
    !     pot2 = real(ueta*area)
    !     pot3 = 0d0
    !     pot4 = 0d0
    !     call dc3d0(alpha, x, y, z, depth, dip, pot1, pot2, pot3, pot4, ux, &
    !                uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz, &
    !                iret)
    !     ! inverse rotation
    !     ux_rot = sin(strike_rad)*ux - cos(strike_rad)*uy
    !     uy_rot = cos(strike_rad)*ux + sin(strike_rad)*uy
    !     ! adjust the unit of variables
    !     ! uret = {pow(10, 2) * ux_rot, pow(10, 2) * uy_rot, pow(10, 2) * uz}
    !     uret(1) = 100d0*ux_rot
    !     uret(2) = 100d0*uy_rot
    !     uret(3) = 100d0*uz
    ! end subroutine

    subroutine call_dc3d0(xsource, ysource, zsource, xobs, yobs, uxi, ueta, &
                          dip, area, strike, uret)
        implicit none
        double precision, intent(in) :: xsource, ysource, zsource, xobs, yobs, &
            uxi, ueta, area, strike
        double precision, intent(inout) :: uret(:)
        double precision, intent(in) :: dip
        double precision :: M_PI, strike_rad, ux_rot, uy_rot
        double precision :: x, y, alpha, z, depth, pot1, pot2, pot3, pot4, &
            ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz
        integer :: iret
        M_PI = 2d0*asin(1d0)
        strike_rad = strike/180d0*M_PI
        ! in DC3D, x axis is the strike direction
        ! rotation
        x = sin(strike_rad)*(xobs - xsource) + cos(strike_rad)*(yobs - ysource)
        y = -cos(strike_rad)*(xobs - xsource) + sin(strike_rad)*(yobs - ysource)
        alpha = (2d0/3d0)
        z = 0d0
        depth = -(zsource)
        pot1 = (uxi*area)
        pot2 = (ueta*area)
        pot3 = 0d0
        pot4 = 0d0
        call dc3d0(alpha, x, y, z, depth, dip, pot1, pot2, pot3, pot4, ux, &
                   uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz, &
                   iret)
        ! inverse rotation
        ux_rot = sin(strike_rad)*ux - cos(strike_rad)*uy
        uy_rot = cos(strike_rad)*ux + sin(strike_rad)*uy
        ! adjust the unit of variables
        ! uret = {pow(10, 2) * ux_rot, pow(10, 2) * uy_rot, pow(10, 2) * uz}
        uret(1) = 100d0*ux_rot
        uret(2) = 100d0*uy_rot
        uret(3) = 100d0*uz
    end subroutine

    subroutine calc_responce(cny_fault, coor_fault, slip_dist, xobs, yobs, xf, &
                             yf, zf, strike, dip, target_id_val, target_id_size, &
                             node_id_in_patch, xinode, etanode, uxinode, uetanode, &
                             r1vec, r2vec, nvec, uobs, uret)
        implicit none
        double precision, intent(in) :: coor_fault(:, :), slip_dist(:, :)
        double precision, intent(in) :: xf, yf, zf, strike, dip, xobs, yobs
        double precision, intent(inout) :: xinode(:), etanode(:), uxinode(:), &
            uetanode(:), r1vec(:), r2vec(:), nvec(:), uobs(:), uret(:)
        integer, intent(in) :: cny_fault(:, :), target_id_val(:)
        integer, intent(in) :: target_id_size
        integer, intent(inout) :: node_id_in_patch(:)
        integer :: itarget, ie, inode, jnode, idim
        integer ::  ir1, ir2
        double precision ::  r1, r2, dip_rad, strike_rad, area, xi, eta, uxi, ueta, &
            xsource, ysource, zsource, M_PI
        M_PI = 2d0*asin(1d0)
        dip_rad = dip/180d0*M_PI
        strike_rad = strike/180d0*M_PI
        ! r1vec(1) = -0.5d0 ! -1d0/sqrt(3d0)
        ! r1vec(2) = 0.5d0 ! 1d0/sqrt(3d0)
        ! r2vec(1) = -0.5d0 ! -1d0/sqrt(3d0)
        ! r2vec(2) = 0.5d0 ! 1d0/sqrt(3d0)
        r1vec(1) = -1d0/sqrt(3d0)
        r1vec(2) = 1d0/sqrt(3d0)
        r2vec(1) = -1d0/sqrt(3d0)
        r2vec(2) = 1d0/sqrt(3d0)
        ! initialize
        do idim = 1, 3
            uobs(idim) = 0d0
        end do

        !   loop for target patchs
        !   (patchs containing the node)
        do itarget = 1, target_id_size
            ie = target_id_val(itarget)
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

            ! surface area of the patch
            area = (xinode(2) - xinode(1))*(etanode(3) - etanode(2))
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
                    do inode = 1, 4
                        xi = xi + nvec(inode)*xinode(inode)
                        eta = eta + nvec(inode)*etanode(inode)
                        uxi = uxi + nvec(inode)*uxinode(inode)
                        ueta = ueta + nvec(inode)*uetanode(inode)
                    end do

                    !   location of the point source specified by
                    !   global coordinates(x, y, z)
                    xsource = xf - eta*cos(dip_rad)*cos(strike_rad) + &
                              xi*sin(strike_rad)
                    ysource = yf + eta*cos(dip_rad)*sin(strike_rad) + &
                              xi*cos(strike_rad)
                    zsource = zf + eta*sin(dip_rad)
                    ! calculate displacement by Okada model
                    call call_dc3d0(xsource, ysource, zsource, xobs, yobs, uxi, ueta, &
                                    dip, area/4d0, strike, uret)
                    ! if (zsource > 0d0) then
                    !     print *, "zsource is positive:", zsource
                    !     print *, "uret: ", uret
                    ! end if
                    !   add contribution from the point source
                    do idim = 1, 3
                        uobs(idim) = uobs(idim) + uret(idim)
                    end do
                end do
            end do
        end do
    end subroutine calc_responce

    subroutine calc_responce_dist(obs_points, cny_fault, coor_fault, slip_dist, &
                                  xf, yf, zf, strike, dip, target_id_val, &
                                  target_id_size, ngnss, node_id_in_patch, &
                                  xinode, etanode, uxinode, uetanode, r1vec, r2vec, &
                                  nvec, response_dist, uobs, uret, nsar_total)
        implicit none
        double precision, intent(in) :: obs_points(:, :), coor_fault(:, :), slip_dist(:, :)
        double precision, intent(in) :: xf, yf, zf, strike, dip
        double precision, intent(inout) :: xinode(:), etanode(:), uxinode(:), &
            uetanode(:), r1vec(:), r2vec(:), nvec(:), response_dist(:, :), uobs(:), uret(:)
        integer, intent(in) :: cny_fault(:, :), target_id_val(:)
        integer, intent(inout) ::  node_id_in_patch(:)
        integer, intent(in) :: ngnss, target_id_size, nsar_total
        integer ::  iobs, idim
        double precision :: xobs, yobs
        ! for SAR observation points
        do iobs = 1, nsar_total
            xobs = obs_points(1, iobs)
            yobs = obs_points(2, iobs)
            ! calculate displacement(ux, uy, uz) at single obsevation point
            call calc_responce(cny_fault, coor_fault, slip_dist, xobs, yobs, xf, &
                               yf, zf, strike, dip, target_id_val, target_id_size, &
                               node_id_in_patch, xinode, etanode, uxinode, uetanode, &
                               r1vec, r2vec, nvec, uobs, uret)
            do idim = 1, 3
                response_dist(idim, iobs) = uobs(idim)
            end do
        end do
        ! for GNSS observation points
        do iobs = 1, ngnss
            xobs = obs_points(1, nsar_total + 1 + 3*(iobs - 1))
            yobs = obs_points(2, nsar_total + 1 + 3*(iobs - 1))
            ! calculate displacement (ux, uy, uz) at single obsevation point
            call calc_responce(cny_fault, coor_fault, slip_dist, xobs, yobs, xf, &
                               yf, zf, strike, dip, target_id_val, target_id_size, &
                               node_id_in_patch, xinode, etanode, uxinode, uetanode, &
                               r1vec, r2vec, nvec, uobs, uret)
            ! copy for three components
            do idim = 1, 3
                response_dist(idim, nsar_total + 3*(iobs - 1) + 1) = uobs(idim)
                response_dist(idim, nsar_total + 3*(iobs - 1) + 2) = uobs(idim)
                response_dist(idim, nsar_total + 3*(iobs - 1) + 3) = uobs(idim)
            end do
        end do
    end subroutine calc_responce_dist

    subroutine calc_greens_func(theta, nplane, nxi_ls, neta_ls, gmat, slip_dist, cny_fault, coor_fault, obs_points, &
                                obs_unitvec, node_to_elem_val, node_to_elem_size, &
                                id_dof, ngnss, nobs, nnode_total, ndof_total, ndof_index, target_id_val, &
                                node_id_in_patch, xinode, etanode, uxinode, uetanode, &
                                r1vec, r2vec, nvec, response_dist, uobs, uret, nsar_total, npath, &
                                nsar_index, gmat_arr, xmin, xmax, zmin, fix_xbend, xbend)
        implicit none
        double precision, intent(inout) :: gmat(:, :), slip_dist(:, :), &
            xinode(:), etanode(:), uxinode(:), uetanode(:), &
            r1vec(:), r2vec(:), nvec(:), response_dist(:, :), uobs(:), uret(:)
        integer, intent(inout) :: target_id_val(:), node_id_in_patch(:)
        integer, intent(in) :: cny_fault(:, :), node_to_elem_val(:, :), &
                               node_to_elem_size(:), id_dof(:), nsar_total
        double precision, intent(in) :: theta(:), coor_fault(:, :), obs_points(:, :), &
            obs_unitvec(:, :), xbend(:)
        integer, intent(in) ::  nplane, nxi_ls(:), neta_ls(:), ngnss, nobs, &
                               nnode_total, ndof_total, ndof_index(:), npath, &
                               nsar_index(:)
        logical, intent(in) :: fix_xbend
        type(mat), intent(inout) :: gmat_arr(:)
        integer :: iplane, idof, inode, idirection, itarget, iobs, idim, i, j, ipath, k
        integer :: target_id_size, nxi, neta
        double precision :: xf, yf, zf, strike, dip, lxi, leta, dip_rad, strike_rad, &
            xmax, xmin, zmin, ymin, ymax, xl, xr, yl, yr
        double precision :: pi

        pi = 4d0*atan(1d0)

        ! initialize gmat
        do j = 1, 2*ndof_total
            do i = 1, nobs
                gmat(i, j) = 0d0
            end do
        end do
        ! initialize slip distribution
        do inode = 1, nnode_total
            do idirection = 1, 2
                slip_dist(idirection, inode) = 0d0
            end do
        end do

        ymin = theta(1)
        ymax = theta(2)
        do iplane = 1, nplane
            call get_geometry(iplane, nplane, theta, nxi, nxi_ls, neta, neta_ls, &
                              xl, yl, xr, yr, xmin, ymin, xmax, ymax, dip, dip_rad, &
                              lxi, leta, zmin, fix_xbend, xbend)

            strike_rad = atan2(xr - xl, yr - yl)
            strike = strike_rad/pi*180d0

            ! xf, yf, zf is coordinate of fault center
            xf = (xl + xr)/2d0 - (zmin*cos(strike_rad))/(2d0*tan(dip_rad))
            yf = (yl + yr)/2d0 + (zmin*sin(strike_rad))/(2d0*tan(dip_rad))
            zf = zmin/2d0
            print *, "---------------"
            print *, xf, yf, zf, lxi, leta
            print *, "---------------"

            ! loop for each degree of freedom of slip
            do idof = ndof_index(iplane), ndof_index(iplane + 1) - 1
                inode = id_dof(idof)
                target_id_size = node_to_elem_size(inode)
                do itarget = 1, target_id_size
                    target_id_val(itarget) = node_to_elem_val(itarget, inode)
                end do
                do idirection = 1, 2
                    ! slip distribution with single unit slip
                    call gen_unit_slip(inode, idirection, slip_dist)
                    ! calculate displacement (x, y, z components) at all the
                    ! observation points
                    call calc_responce_dist(obs_points, cny_fault, coor_fault, slip_dist, &
                                            xf, yf, zf, strike, dip, target_id_val, &
                                            target_id_size, ngnss, node_id_in_patch, &
                                            xinode, etanode, uxinode, uetanode, r1vec, r2vec, &
                                            nvec, response_dist, uobs, uret, nsar_total)

                    ! inner product (displacement * LOS unitvec)
                    do iobs = 1, nobs
                        do idim = 1, 3
                            gmat(iobs, 2*(idof - 1) + idirection) = &
                                gmat(iobs, 2*(idof - 1) + idirection) + &
                                response_dist(idim, iobs)*obs_unitvec(idim, iobs)
                        end do
                    end do
                    ! slip distribution with single unit slip
                    call del_unit_slip(inode, idirection, slip_dist)
                end do
            end do
        end do

        do ipath = 1, npath
            k = nsar_index(ipath)
            do iobs = nsar_index(ipath), nsar_index(ipath + 1) - 1
                do idim = 1, 2*ndof_total
                    gmat_arr(ipath)%body(iobs - k + 1, idim) = gmat(iobs, idim)
                end do
            end do
            ! print *, "ipath: ", ipath
            ! do iobs = nsar_index(ipath), nsar_index(ipath + 1) - 1
            !     do idim = 1, 2*ndof_total
            !         write(*, "(e15.5)", advance="no") gmat_arr(ipath)%body(iobs - nsar_index(ipath) + 1, idim)
            !     end do
            !     write(*, *)
            ! end do
        end do
    end subroutine
end module gfunc
