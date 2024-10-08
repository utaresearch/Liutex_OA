module liutex_mod
    !!-------------------------------------------------------------------------
    !! Module that contains subroutines and functions to aid in the 
    !! research of Liutex.
    !!
    !! By: Oscar Alvarez
    !! email: oscar.alvarez@uta.edu
    !!-------------------------------------------------------------------------
    implicit none

    !! Interfaces


    !! Parameters

    contains

    !! Subroutines
    subroutine liutex(a, vor, liutex_vec, liutex_mag)
        implicit none

        real(8), dimension(3,3), intent(in) :: a            ! Velocity gradient tensor
        real(8), dimension(3), intent(in) :: vor            ! Vorticity
        real(8), dimension(3), intent(out) :: liutex_vec    ! Liutex direction vector 
        real(8), intent(out) :: liutex_mag                  ! Liutex magnitude
      
        real(8), dimension(3,3) :: tt
        real(8) :: aa, bb, cc, delta, rr, aaaa, bbbb, qq, delta1, delta2, delta3, temp
        real(8) :: eig3r
        complex(8) :: eig1c, eig2c 

        ! Cubic Formula
        ! Reference: Numerical Recipes in FORTRAN 77, Second Edition
        ! 5.6 Quadratic and Cubic Equations
        ! Page 179
        !---------------------------------------------------------------------

        ! cubic equation
        ! x**3 + aa * x**2 + bb * x + cc = 0

        ! coefficients of characteristic equation

        aa = -(a(1,1) + a(2,2) + a(3,3))

        tt = matmul(a,a)

        bb = -0.5d0*(tt(1,1) + tt(2,2) + tt(3,3) - (a(1,1) + a(2,2) + a(3,3))**2)

        cc = -( a(1,1) * (a(2,2)*a(3,3)-a(2,3)*a(3,2))                &
                - a(1,2) * (a(2,1)*a(3,3)-a(2,3)*a(3,1))              &
                + a(1,3) * (a(2,1)*a(3,2)-a(2,2)*a(3,1)) )

        ! discriminant of characteristic equation
        delta = 18.d0*aa*bb*cc - 4.d0*aa**3 * cc + aa**2 * bb**2 - 4.d0*bb**3 - 27.d0*cc**2

        qq = (aa**2.d0 - 3.d0*bb) / 9.d0
        rr = (2.d0*aa**3 - 9.d0*aa*bb + 27.d0*cc) / 54.d0

        ! alleviate round error
        delta = -delta / 108.d0
        liutex_vec = 0.d0
        liutex_mag = 0.d0

        if (delta > 0.d0) then ! one real root and two complex conjugate roots

            aaaa = -sign(1.d0, rr) * (abs(rr)+sqrt(delta))**(1.d0/3.d0)

            if (aaaa == 0.d0) then
                bbbb = 0.d0
            else
                bbbb = qq / aaaa
            end if

            eig1c = cmplx(-0.5*(aaaa+bbbb)-aa/3.0, 0.5*sqrt(3.0)*(aaaa-bbbb), kind=8)
            eig2c = cmplx(real(eig1c,kind=8), -aimag(eig1c), kind=8) 
            eig3r = aaaa + bbbb - aa/3.d0

            ! real right eigenvector

            delta1 = (a(1,1) - eig3r) * (a(2,2) - eig3r) - a(2,1)*a(1,2)
            delta2 = (a(2,2) - eig3r) * (a(3,3) - eig3r) - a(2,3)*a(3,2)
            delta3 = (a(1,1) - eig3r) * (a(3,3) - eig3r) - a(1,3)*a(3,1)

            if (delta1 == 0.d0 .and. delta2 == 0.d0 .and. delta3 == 0.d0) then
                write(*,*) 'ERROR: delta1 = delta2 = delta3 = 0.0'
                write(*,*) a(1,1)-eig3r,  a(1,2),       a(1,3)
                write(*,*) a(2,1),        a(2,2)-eig3r, a(2,3)
                write(*,*) a(3,1),        a(3,2),       a(3,3)-eig3r
                stop
            end if

            if (abs(delta1) >= abs(delta2) .and. abs(delta1) >= abs(delta3)) then

                liutex_vec(1) = (-(a(2,2)-eig3r)*a(1,3) +         a(1,2)*a(2,3))/delta1
                liutex_vec(2) = (         a(2,1)*a(1,3) - (a(1,1)-eig3r)*a(2,3))/delta1
                liutex_vec(3) = 1.d0

            else if (abs(delta2) >= abs(delta1) .and. abs(delta2) >= abs(delta3)) then

                liutex_vec(1) = 1.d0
                liutex_vec(2) = (-(a(3,3)-eig3r)*a(2,1) +         a(2,3)*a(3,1))/delta2
                liutex_vec(3) = (         a(3,2)*a(2,1) - (a(2,2)-eig3r)*a(3,1))/delta2

            else if (abs(delta3) >= abs(delta1) .and. abs(delta3) >= abs(delta2)) then

                liutex_vec(1) = (-(a(3,3)-eig3r)*a(1,2) +         a(1,3)*a(3,2))/delta3
                liutex_vec(2) = 1.d0
                liutex_vec(3) = (         a(3,1)*a(1,2) - (a(1,1)-eig3r)*a(3,2))/delta3

            else
                write(*,*) 'ERROR: '
                write(*,*) delta1, delta2, delta3
                stop
            end if

            temp = sqrt(liutex_vec(1)**2 + liutex_vec(2)**2 + liutex_vec(3)**2)

            liutex_vec(1) = liutex_vec(1) / temp
            liutex_vec(2) = liutex_vec(2) / temp
            liutex_vec(3) = liutex_vec(3) / temp

            temp = dot_product(vor, liutex_vec)
            
            liutex_vec(1) = sign(1.d0,temp) * liutex_vec(1)
            liutex_vec(2) = sign(1.d0,temp) * liutex_vec(2)
            liutex_vec(3) = sign(1.d0,temp) * liutex_vec(3)

            liutex_mag = (temp - sqrt(temp**2 - 4.d0*aimag(eig2c)*aimag(eig2c)))
        end if
      
    end subroutine liutex

    !! Functions

    function vorticity(nabla_v) result(vor)
        !!! Calculate vorticity vector from velocity gradient tensor
        implicit none
        real(8), dimension(3,3), intent(in) :: nabla_v
        real(8), dimension(3) :: vor

        vor = 0.d0
        vor(1) = nabla_v(3,2) - nabla_v(2,3)
        vor(2) = nabla_v(1,3) - nabla_v(3,1)
        vor(3) = nabla_v(2,1) - nabla_v(1,2)
    end function vorticity


    function velocity_gradient(u, v, w, x, y, z, imax, jmax, kmax) result(nabla_v)
        !!! Calculates the velocity gradient tensor using the u, v, and w velocities.
        implicit none
        
        integer, intent(in) :: imax, jmax, kmax
        real(8), dimension(imax,jmax,kmax), intent(in) :: u, v, w, x, y, z
        real(8), dimension(imax,jmax,kmax,3,3) :: nabla_v

        real(8) :: u_xi, u_eta, u_zeta
        real(8) :: v_xi, v_eta, v_zeta
        real(8) :: w_xi, w_eta, w_zeta
        real(8) :: x_xi, x_eta, x_zeta
        real(8) :: y_xi, y_eta, y_zeta
        real(8) :: z_xi, z_eta, z_zeta
        real(8) :: xi_x, xi_y, xi_z
        real(8) :: eta_x, eta_y, eta_z
        real(8) :: zeta_x, zeta_y, zeta_z
        real(8) :: det
        integer :: i, j, k

        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax
        
                    ! Using Finite difference scheme to calculate partial derivatives

                    if(i == 1) then
                        ! forward difference
                        u_xi = u(2, j, k) - u(1, j, k)
                        v_xi = v(2, j, k) - v(1, j, k)
                        w_xi = w(2, j, k) - w(1, j, k)
            
                        x_xi = x(2, j, k) - x(1, j, k)
                        y_xi = y(2, j, k) - y(1, j, k)
                        z_xi = z(2, j, k) - z(1, j, k)
                    else if(i == imax) then
                        ! backward difference
                        u_xi = u(imax, j, k) - u(imax-1, j, k)
                        v_xi = v(imax, j, k) - v(imax-1, j, k)
                        w_xi = w(imax, j, k) - w(imax-1, j, k)
            
                        x_xi = x(imax, j, k) - x(imax-1, j, k)
                        y_xi = y(imax, j, k) - y(imax-1, j, k)
                        z_xi = z(imax, j, k) - z(imax-1, j, k)
                    else
                        ! central difference
                        u_xi = 0.5*(u(i+1, j, k) - u(i-1, j, k))
                        v_xi = 0.5*(v(i+1, j, k) - v(i-1, j, k))
                        w_xi = 0.5*(w(i+1, j, k) - w(i-1, j, k))
            
                        x_xi = 0.5*(x(i+1, j, k) - x(i-1, j, k))
                        y_xi = 0.5*(y(i+1, j, k) - y(i-1, j, k))
                        z_xi = 0.5*(z(i+1, j, k) - z(i-1, j, k))
                    end if
        
                    if(j == 1) then
                        ! forward difference
                        u_eta = u(i, 2, k) - u(i, 1, k)
                        v_eta = v(i, 2, k) - v(i, 1, k)
                        w_eta = w(i, 2, k) - w(i, 1, k)
            
                        x_eta = x(i, 2, k) - x(i, 1, k)
                        y_eta = y(i, 2, k) - y(i, 1, k)
                        z_eta = z(i, 2, k) - z(i, 1, k)
                    else if(j == jmax) then
                        ! backward difference
                        u_eta = u(i, jmax, k) - u(i, jmax-1, k)
                        v_eta = v(i, jmax, k) - v(i, jmax-1, k)
                        w_eta = w(i, jmax, k) - w(i, jmax-1, k)
            
                        x_eta = x(i, jmax, k) - x(i, jmax-1, k)
                        y_eta = y(i, jmax, k) - y(i, jmax-1, k)
                        z_eta = z(i, jmax, k) - z(i, jmax-1, k)    
                    else
                        ! central difference
                        u_eta = 0.5*(u(i, j+1, k) - u(i, j-1, k))
                        v_eta = 0.5*(v(i, j+1, k) - v(i, j-1, k))
                        w_eta = 0.5*(w(i, j+1, k) - w(i, j-1, k))
            
                        x_eta = 0.5*(x(i, j+1, k) - x(i, j-1, k))
                        y_eta = 0.5*(y(i, j+1, k) - y(i, j-1, k))
                        z_eta = 0.5*(z(i, j+1, k) - z(i, j-1, k))    
                    end if
        
                    if(k == 1) then    
                        ! forward difference
                        u_zeta = u(i, j, 2) - u(i, j, 1)
                        v_zeta = v(i, j, 2) - v(i, j, 1)
                        w_zeta = w(i, j, 2) - w(i, j, 1)
            
                        x_zeta = x(i, j, 2) - x(i, j, 1)
                        y_zeta = y(i, j, 2) - y(i, j, 1)
                        z_zeta = z(i, j, 2) - z(i, j, 1)
                    else if(k == kmax) then
                        ! backward difference
                        u_zeta = u(i, j, kmax) - u(i, j, kmax-1)
                        v_zeta = v(i, j, kmax) - v(i, j, kmax-1)
                        w_zeta = w(i, j, kmax) - w(i, j, kmax-1)
            
                        x_zeta = x(i, j, kmax) - x(i, j, kmax-1)
                        y_zeta = y(i, j, kmax) - y(i, j, kmax-1)
                        z_zeta = z(i, j, kmax) - z(i, j, kmax-1)
                    else
                        ! central difference
                        u_zeta = 0.5*(u(i, j, k+1) - u(i, j, k-1))
                        v_zeta = 0.5*(v(i, j, k+1) - v(i, j, k-1))
                        w_zeta = 0.5*(w(i, j, k+1) - w(i, j, k-1))
            
                        x_zeta = 0.5*(x(i, j, k+1) - x(i, j, k-1))
                        y_zeta = 0.5*(y(i, j, k+1) - y(i, j, k-1))
                        z_zeta = 0.5*(z(i, j, k+1) - z(i, j, k-1))
                    end if
        
                    ! determinant of Jacobian
                    det =   x_xi * (y_eta*z_zeta-y_zeta*z_eta)                               &
                            - x_eta * (y_xi*z_zeta-y_zeta*z_xi)                              &
                            + x_zeta * (y_xi*z_eta-y_eta*z_xi)
        
                    det = 1.0 / det
        
                    xi_x = det*(y_eta*z_zeta - y_zeta*z_eta)
                    xi_y = det*(x_zeta*z_eta - x_eta*z_zeta)
                    xi_z = det*(x_eta*y_zeta - x_zeta*y_eta)
        
                    eta_x = det*(y_zeta*z_xi - y_xi*z_zeta)
                    eta_y = det*(x_xi*z_zeta - x_zeta*z_xi)
                    eta_z = det*(x_zeta*y_xi - x_xi*y_zeta)
        
                    zeta_x = det*(y_xi*z_eta - y_eta*z_xi)
                    zeta_y = det*(x_eta*z_xi - x_xi*z_eta)
                    zeta_z = det*(x_xi*y_eta - x_eta*y_xi)

                    ! Assembling the velocity gradient tensor
                    
                    nabla_v(i, j, k, 1, 1) = u_xi*xi_x + u_eta*eta_x + u_zeta*zeta_x    !dudx
                    nabla_v(i, j, k, 1, 2) = u_xi*xi_y + u_eta*eta_y + u_zeta*zeta_y    !dudy
                    nabla_v(i, j, k, 1, 3) = u_xi*xi_z + u_eta*eta_z + u_zeta*zeta_z    !dudz
        
                    nabla_v(i, j, k, 2, 1) = v_xi*xi_x + v_eta*eta_x + v_zeta*zeta_x    !dvdx
                    nabla_v(i, j, k, 2, 2) = v_xi*xi_y + v_eta*eta_y + v_zeta*zeta_y    !dvdy
                    nabla_v(i, j, k, 2, 3) = v_xi*xi_z + v_eta*eta_z + v_zeta*zeta_z    !dvdz
        
                    nabla_v(i, j, k, 3, 1) = w_xi*xi_x + w_eta*eta_x + w_zeta*zeta_x    !dwdx
                    nabla_v(i, j, k, 3, 2) = w_xi*xi_y + w_eta*eta_y + w_zeta*zeta_y    !dwdy
                    nabla_v(i, j, k, 3, 3) = w_xi*xi_z + w_eta*eta_z + w_zeta*zeta_z    !dwdz
                
                end do
            end do
        end do

    end function velocity_gradient


    function interpolation(f_x, f_y, f_z, x, y, z, i, j, k, imax, jmax, kmax) result(big_f)
        !!!!!! UNFINISHED !!!!!!!!!!!
        !! Integrates a vector valued function f = (f_x, f_y, f_z) by the method described
        !! in Fudan paper.
        implicit none
        integer, intent(in) :: i, j, k, imax, jmax, kmax
        real(8), dimension(imax,jmax,kmax), intent(in) :: f_x, f_y, f_z, x, y, z
        real(8), dimension(imax,jmax,kmax,3) :: big_f

        real(8), dimension(8) :: alpha
        real(8), dimension(8) :: x_s, y_s, z_s, fx_s, fy_s, fz_s
        real(8), dimension(3) :: ds
        real(8), dimension(3) :: interpolated_sum
        real(8) :: x_i, y_i, z_i
        real(8) :: fx_i, fy_i, fz_i
        logical :: x_boundary_stop, y_boundary_stop, z_boundary_stop, boundary_stop
        logical :: function_value_stop
        integer :: max_n
        integer :: n, s

        !! Max number of iterations/generated seed points
        max_n = 100

        !! Distance to be incremented when interpolating values.
        ds(1) = (x(i+1,j+1,k+1) - x(i,j,k)) / 10.d0
        ds(2) = (y(i+1,j+1,k+1) - y(i,j,k)) / 10.d0
        ds(3) = (z(i+1,j+1,k+1) - z(i,j,k)) / 10.d0

        !! Building Volume control unit. Each index is a vertex of the control unit.
        x_s(1) = x(i,j,k)
        x_s(2) = x(i+1,j,k)
        x_s(3) = x(i+1,j+1,k)
        x_s(4) = x(i,j+1,k)
        x_s(5) = x(i,j,k+1)
        x_s(6) = x(i+1,j,k+1)
        x_s(7) = x(i+1,j+1,k+1)
        x_s(8) = x(i,j+1,k+1)

        y_s(1) = y(i,j,k)
        y_s(2) = y(i+1,j,k)
        y_s(3) = y(i+1,j+1,k)
        y_s(4) = y(i,j+1,k)
        y_s(5) = y(i,j,k+1)
        y_s(6) = y(i+1,j,k+1)
        y_s(7) = y(i+1,j+1,k+1)
        y_s(8) = y(i,j+1,k+1)

        z_s(1) = z(i,j,k)
        z_s(2) = z(i+1,j,k)
        z_s(3) = z(i+1,j+1,k)
        z_s(4) = z(i,j+1,k)
        z_s(5) = z(i,j,k+1)
        z_s(6) = z(i+1,j,k+1)
        z_s(7) = z(i+1,j+1,k+1)
        z_s(8) = z(i,j+1,k+1)

        fx_s(1) = f_x(i,j,k)
        fx_s(2) = f_x(i+1,j,k)
        fx_s(3) = f_x(i+1,j+1,k)
        fx_s(4) = f_x(i,j+1,k)
        fx_s(5) = f_x(i,j,k+1)
        fx_s(6) = f_x(i+1,j,k+1)
        fx_s(7) = f_x(i+1,j+1,k+1)
        fx_s(8) = f_x(i,j+1,k+1)

        fy_s(1) = f_y(i,j,k)
        fy_s(2) = f_y(i+1,j,k)
        fy_s(3) = f_y(i+1,j+1,k)
        fy_s(4) = f_y(i,j+1,k)
        fy_s(5) = f_y(i,j,k+1)
        fy_s(6) = f_y(i+1,j,k+1)
        fy_s(7) = f_y(i+1,j+1,k+1)
        fy_s(8) = f_y(i,j+1,k+1)

        fz_s(1) = f_z(i,j,k)
        fz_s(2) = f_z(i+1,j,k)
        fz_s(3) = f_z(i+1,j+1,k)
        fz_s(4) = f_z(i,j+1,k)
        fz_s(5) = f_z(i,j,k+1)
        fz_s(6) = f_z(i+1,j,k+1)
        fz_s(7) = f_z(i+1,j+1,k+1)
        fz_s(8) = f_z(i,j+1,k+1)

        !! Initializing values (current position/value)
        x_i = x(i,j,k)
        y_i = y(i,j,k)
        z_i = z(i,j,k)

        fx_i = f_x(i,j,k)
        fy_i = f_y(i,j,k)
        fz_i = f_y(i,j,k)

        interpolated_sum = 0.d0

        do n = 1, max_n

            !! Stop if we are outside of the calculation domain (control unit).
            x_boundary_stop = (x_i < x(i,j,k)) .or. (x_i > x(i,j,k))
            y_boundary_stop = (y_i < y(i,j,k)) .or. (y_i > y(i,j,k))
            z_boundary_stop = (z_i < z(i,j,k)) .or. (z_i > z(i,j,k))
            boundary_stop = x_boundary_stop .or. y_boundary_stop .or. z_boundary_stop

            !! Stop if function value is = 0.
            function_value_stop = (norm2((/fx_i, fy_i, fz_i/)) == 0.d0)

            if (boundary_stop .or. function_value_stop) then
                exit
            end if

            !! Interpolate function value
            do s = 1, 8
                alpha(s) = (x_s(7) - x_i) * (y_s(7) - y_i) * (z_s(7) - z_i)
                alpha(s) = alpha(s) / ( (x_s(7) - x_s(s))*(y_s(7)-y_s(s))*(z_s(7)-z_s(s)) )

                interpolated_sum(1) = interpolated_sum(1) + alpha(s)*fx_s(s)
                interpolated_sum(2) = interpolated_sum(2) + alpha(s)*fy_s(s)
                interpolated_sum(3) = interpolated_sum(3) + alpha(s)*fz_s(s)
            end do
            
            !! Interpolated values of function.
            fx_i = interpolated_sum(1)
            fy_i = interpolated_sum(2)
            fz_i = interpolated_sum(3)

            !! Position/location of interpolated function value.
            x_i = x_i + ds(1)*fx_i
            y_i = y_i + ds(2)*fy_i
            z_i = z_i + ds(3)*fz_i
            
        end do
            


    end function interpolation


    function hessian_mat(f, x, y, z, imax, jmax, kmax) result(h_mat)
        !! Creates the Hessian Matrix for each node in 3d
        implicit none
        integer, intent(in) :: imax, jmax, kmax
        real(8), dimension(imax,jmax,kmax), intent(in) :: f, x, y, z
        real(8), dimension(imax,jmax,kmax,3,3) :: h_mat

        real(8), dimension(imax,jmax,kmax,3) :: df
        real(8), dimension(imax,jmax,kmax) :: f_x, f_y, f_z
        real(8), dimension(imax,jmax,kmax) :: xi_x, xi_y, xi_z
        real(8), dimension(imax,jmax,kmax) :: eta_x, eta_y, eta_z
        real(8), dimension(imax,jmax,kmax) :: zeta_x, zeta_y, zeta_z

        real(8) :: fx_xi, fx_eta, fx_zeta
        real(8) :: fy_xi, fy_eta, fy_zeta
        real(8) :: fz_xi, fz_eta, fz_zeta
        real(8) :: x_xi, x_eta, x_zeta
        real(8) :: y_xi, y_eta, y_zeta
        real(8) :: z_xi, z_eta, z_zeta
        real(8) :: dfx_dx, dfx_dy, dfx_dz
        real(8) :: dfy_dx, dfy_dy, dfy_dz
        real(8) :: dfz_dx, dfz_dy, dfz_dz
        
        real(8) :: det
        integer :: i, j, k

        df = gradient(f, x, y, z, imax, jmax, kmax)

        f_x = df(:,:,:,1)
        f_y = df(:,:,:,2)
        f_z = df(:,:,:,3)

        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax

                    !! Finite-Difference in parametric space.
                    if(i == 1) then
                        fx_xi = f_x(2, j, k) - f_x(1, j, k)
                        fy_xi = f_y(2, j, k) - f_y(1, j, k)
                        fz_xi = f_z(2, j, k) - f_z(1, j, k)
                        
                        x_xi = x(2, j, k) - x(1, j, k)
                        y_xi = y(2, j, k) - y(1, j, k)
                        z_xi = z(2, j, k) - z(1, j, k)
                    else if(i == imax) then
                        fx_xi = f_x(imax, j, k) - f_x(imax-1, j, k)
                        fy_xi = f_y(imax, j, k) - f_y(imax-1, j, k)
                        fz_xi = f_z(imax, j, k) - f_z(imax-1, j, k)
                        
                        x_xi = x(imax, j, k)-x(imax-1, j, k)
                        y_xi = y(imax, j, k)-y(imax-1, j, k)
                        z_xi = z(imax, j, k)-z(imax-1, j, k)
                    else
                        fx_xi = 0.5d0*(f_x(i+1, j, k) - f_x(i-1, j, k))
                        fy_xi = 0.5d0*(f_y(i+1, j, k) - f_y(i-1, j, k))
                        fz_xi = 0.5d0*(f_z(i+1, j, k) - f_z(i-1, j, k))

                        x_xi = 0.5d0*(x(i+1, j, k) - x(i-1, j, k))
                        y_xi = 0.5d0*(y(i+1, j, k) - y(i-1, j, k))
                        z_xi = 0.5d0*(z(i+1, j, k) - z(i-1, j, k))
                    end if

                    if(j == 1) then
                        fx_eta = f_x(i, 2, k) - f_x(i, 1, k)
                        fy_eta = f_y(i, 2, k) - f_y(i, 1, k)
                        fz_eta = f_z(i, 2, k) - f_z(i, 1, k)

                        x_eta = x(i, 2, k) - x(i, 1, k)
                        y_eta = y(i, 2, k) - y(i, 1, k)
                        z_eta = z(i, 2, k) - z(i, 1, k)
                    else if(j == jmax) then
                        fx_eta = f_x(i, jmax, k) - f_x(i, jmax-1, k)
                        fy_eta = f_y(i, jmax, k) - f_y(i, jmax-1, k)
                        fz_eta = f_z(i, jmax, k) - f_z(i, jmax-1, k)

                        x_eta = x(i, jmax, k) - x(i, jmax-1, k)
                        y_eta = y(i, jmax, k) - y(i, jmax-1, k)
                        z_eta = z(i, jmax, k) - z(i, jmax-1, k)
                    else
                        fx_eta = 0.5d0*(f_x(i, j+1, k) - f_x(i, j-1, k))
                        fy_eta = 0.5d0*(f_y(i, j+1, k) - f_y(i, j-1, k))
                        fz_eta = 0.5d0*(f_z(i, j+1, k) - f_z(i, j-1, k))

                        x_eta = 0.5d0*(x(i, j+1, k) - x(i, j-1, k))
                        y_eta = 0.5d0*(y(i, j+1, k) - y(i, j-1, k))
                        z_eta = 0.5d0*(z(i, j+1, k) - z(i, j-1, k))
                    end if

                    if(k == 1) then
                        fx_zeta = f_x(i, j, 2) - f_x(i, j, 1)
                        fy_zeta = f_y(i, j, 2) - f_y(i, j, 1)
                        fz_zeta = f_z(i, j, 2) - f_z(i, j, 1)

                        x_zeta = x(i, j, 2) - x(i, j, 1)
                        y_zeta = y(i, j, 2) - y(i, j, 1)
                        z_zeta = z(i, j, 2) - z(i, j, 1)
                    else if(k == kmax) then
                        fx_zeta = f_x(i, j, kmax) - f_x(i, j, kmax-1)
                        fy_zeta = f_y(i, j, kmax) - f_y(i, j, kmax-1)
                        fz_zeta = f_z(i, j, kmax) - f_z(i, j, kmax-1)

                        x_zeta = x(i, j, kmax) - x(i, j, kmax-1)
                        y_zeta = y(i, j, kmax) - y(i, j, kmax-1)
                        z_zeta = z(i, j, kmax) - z(i, j, kmax-1)
                    else
                        fx_zeta = 0.5d0*(f_x(i, j, k+1) - f_x(i, j, k-1))
                        fy_zeta = 0.5d0*(f_y(i, j, k+1) - f_y(i, j, k-1))
                        fz_zeta = 0.5d0*(f_z(i, j, k+1) - f_z(i, j, k-1))

                        x_zeta = 0.5d0*(x(i, j, k+1) - x(i, j, k-1))
                        y_zeta = 0.5d0*(y(i, j, k+1) - y(i, j, k-1))
                        z_zeta = 0.5d0*(z(i, j, k+1) - z(i, j, k-1))
                    end if

                    !! Jacobian Transformation from parametric to global space.
                    det =   x_xi*(y_eta*z_zeta-y_zeta*z_eta)  &
                            - x_eta*(y_xi*z_zeta-y_zeta*z_xi)   &
                            + x_zeta*(y_xi*z_eta-y_eta*z_xi)

                    det = 1.d0 / det

                    xi_x(i,j,k) = det*(y_eta*z_zeta - y_zeta*z_eta)
                    xi_y(i,j,k) = det*(x_zeta*z_eta - x_eta*z_zeta)
                    xi_z(i,j,k) = det*(x_eta*y_zeta - x_zeta*y_eta)

                    eta_x(i,j,k) = det*(y_zeta*z_xi - y_xi*z_zeta)
                    eta_y(i,j,k) = det*(x_xi*z_zeta - x_zeta*z_xi)
                    eta_z(i,j,k) = det*(x_zeta*y_xi - x_xi*y_zeta)

                    zeta_x(i,j,k) = det*(y_xi*z_eta - y_eta*z_xi)
                    zeta_y(i,j,k) = det*(x_eta*z_xi - x_xi*z_eta)
                    zeta_z(i,j,k) = det*(x_xi*y_eta - x_eta*y_xi)

                    !! Second partial derivatives in global space
                    dfx_dx = fx_xi*xi_x(i,j,k) + fx_eta*eta_x(i,j,k) + fx_zeta*zeta_x(i,j,k)
                    dfx_dy = fx_xi*xi_y(i,j,k) + fx_eta*eta_y(i,j,k) + fx_zeta*zeta_y(i,j,k)
                    dfx_dz = fx_xi*xi_z(i,j,k) + fx_eta*eta_z(i,j,k) + fx_zeta*zeta_z(i,j,k)

                    dfy_dx = fy_xi*xi_x(i,j,k) + fy_eta*eta_x(i,j,k) + fy_zeta*zeta_x(i,j,k)
                    dfy_dy = fy_xi*xi_y(i,j,k) + fy_eta*eta_y(i,j,k) + fy_zeta*zeta_y(i,j,k)
                    dfy_dz = fy_xi*xi_z(i,j,k) + fy_eta*eta_z(i,j,k) + fy_zeta*zeta_z(i,j,k)

                    dfz_dx = fz_xi*xi_x(i,j,k) + fz_eta*eta_x(i,j,k) + fz_zeta*zeta_x(i,j,k)
                    dfz_dy = fz_xi*xi_y(i,j,k) + fz_eta*eta_y(i,j,k) + fz_zeta*zeta_y(i,j,k)
                    dfz_dz = fz_xi*xi_z(i,j,k) + fz_eta*eta_z(i,j,k) + fz_zeta*zeta_z(i,j,k)

                    !! Forming the Hessian matrix
                    h_mat(i,j,k,1,1) = dfx_dx
                    h_mat(i,j,k,1,2) = dfx_dy
                    h_mat(i,j,k,1,3) = dfx_dz
                    h_mat(i,j,k,2,1) = dfy_dx
                    h_mat(i,j,k,2,2) = dfy_dy
                    h_mat(i,j,k,2,3) = dfy_dz
                    h_mat(i,j,k,3,1) = dfz_dx
                    h_mat(i,j,k,3,2) = dfz_dy
                    h_mat(i,j,k,3,3) = dfz_dz
                    
                end do
            end do
        end do

    end function hessian_mat


    function divergence(f_x, f_y, f_z, x, y, z, imax, jmax, kmax) result(div_f)
        !! Finds the divergence of a vector function f = < f_x, f_y, f_z >
        implicit none
        integer, intent(in) :: imax, jmax, kmax
        real(8), dimension(imax,jmax,kmax), intent(in) :: f_x, f_y, f_z, x, y, z
        real(8), dimension(imax,jmax,kmax) :: div_f

        real(8), dimension(imax,jmax,kmax) :: xi_x, xi_y, xi_z
        real(8), dimension(imax,jmax,kmax) :: eta_x, eta_y, eta_z
        real(8), dimension(imax,jmax,kmax) :: zeta_x, zeta_y, zeta_z

        real(8) :: fx_xi, fx_eta, fx_zeta
        real(8) :: fy_xi, fy_eta, fy_zeta
        real(8) :: fz_xi, fz_eta, fz_zeta
        real(8) :: x_xi, x_eta, x_zeta
        real(8) :: y_xi, y_eta, y_zeta
        real(8) :: z_xi, z_eta, z_zeta
        real(8) :: dfx_dx, dfy_dy, dfz_dz
        
        real(8) :: det
        integer :: i, j, k

        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax

                    !! Finite-Difference in parametric space.
                    if(i == 1) then
                        fx_xi = f_x(2, j, k) - f_x(1, j, k)
                        fy_xi = f_y(2, j, k) - f_y(1, j, k)
                        fz_xi = f_z(2, j, k) - f_z(1, j, k)
                        
                        x_xi = x(2, j, k) - x(1, j, k)
                        y_xi = y(2, j, k) - y(1, j, k)
                        z_xi = z(2, j, k) - z(1, j, k)
                    else if(i == imax) then
                        fx_xi = f_x(imax, j, k) - f_x(imax-1, j, k)
                        fy_xi = f_y(imax, j, k) - f_y(imax-1, j, k)
                        fz_xi = f_z(imax, j, k) - f_z(imax-1, j, k)
                        
                        x_xi = x(imax, j, k)-x(imax-1, j, k)
                        y_xi = y(imax, j, k)-y(imax-1, j, k)
                        z_xi = z(imax, j, k)-z(imax-1, j, k)
                    else
                        fx_xi = 0.5d0*(f_x(i+1, j, k) - f_x(i-1, j, k))
                        fy_xi = 0.5d0*(f_y(i+1, j, k) - f_y(i-1, j, k))
                        fz_xi = 0.5d0*(f_z(i+1, j, k) - f_z(i-1, j, k))

                        x_xi = 0.5d0*(x(i+1, j, k) - x(i-1, j, k))
                        y_xi = 0.5d0*(y(i+1, j, k) - y(i-1, j, k))
                        z_xi = 0.5d0*(z(i+1, j, k) - z(i-1, j, k))
                    end if

                    if(j == 1) then
                        fx_eta = f_x(i, 2, k) - f_x(i, 1, k)
                        fy_eta = f_y(i, 2, k) - f_y(i, 1, k)
                        fz_eta = f_z(i, 2, k) - f_z(i, 1, k)

                        x_eta = x(i, 2, k) - x(i, 1, k)
                        y_eta = y(i, 2, k) - y(i, 1, k)
                        z_eta = z(i, 2, k) - z(i, 1, k)
                    else if(j == jmax) then
                        fx_eta = f_x(i, jmax, k) - f_x(i, jmax-1, k)
                        fy_eta = f_y(i, jmax, k) - f_y(i, jmax-1, k)
                        fz_eta = f_z(i, jmax, k) - f_z(i, jmax-1, k)

                        x_eta = x(i, jmax, k) - x(i, jmax-1, k)
                        y_eta = y(i, jmax, k) - y(i, jmax-1, k)
                        z_eta = z(i, jmax, k) - z(i, jmax-1, k)
                    else
                        fx_eta = 0.5d0*(f_x(i, j+1, k) - f_x(i, j-1, k))
                        fy_eta = 0.5d0*(f_y(i, j+1, k) - f_y(i, j-1, k))
                        fz_eta = 0.5d0*(f_z(i, j+1, k) - f_z(i, j-1, k))

                        x_eta = 0.5d0*(x(i, j+1, k) - x(i, j-1, k))
                        y_eta = 0.5d0*(y(i, j+1, k) - y(i, j-1, k))
                        z_eta = 0.5d0*(z(i, j+1, k) - z(i, j-1, k))
                    end if

                    if(k == 1) then
                        fx_zeta = f_x(i, j, 2) - f_x(i, j, 1)
                        fy_zeta = f_y(i, j, 2) - f_y(i, j, 1)
                        fz_zeta = f_z(i, j, 2) - f_z(i, j, 1)

                        x_zeta = x(i, j, 2) - x(i, j, 1)
                        y_zeta = y(i, j, 2) - y(i, j, 1)
                        z_zeta = z(i, j, 2) - z(i, j, 1)
                    else if(k == kmax) then
                        fx_zeta = f_x(i, j, kmax) - f_x(i, j, kmax-1)
                        fy_zeta = f_y(i, j, kmax) - f_y(i, j, kmax-1)
                        fz_zeta = f_z(i, j, kmax) - f_z(i, j, kmax-1)

                        x_zeta = x(i, j, kmax) - x(i, j, kmax-1)
                        y_zeta = y(i, j, kmax) - y(i, j, kmax-1)
                        z_zeta = z(i, j, kmax) - z(i, j, kmax-1)
                    else
                        fx_zeta = 0.5d0*(f_x(i, j, k+1) - f_x(i, j, k-1))
                        fy_zeta = 0.5d0*(f_y(i, j, k+1) - f_y(i, j, k-1))
                        fz_zeta = 0.5d0*(f_z(i, j, k+1) - f_z(i, j, k-1))

                        x_zeta = 0.5d0*(x(i, j, k+1) - x(i, j, k-1))
                        y_zeta = 0.5d0*(y(i, j, k+1) - y(i, j, k-1))
                        z_zeta = 0.5d0*(z(i, j, k+1) - z(i, j, k-1))
                    end if

                    !! Jacobian Transformation from parametric to global space.
                    det =   x_xi*(y_eta*z_zeta-y_zeta*z_eta)  &
                            - x_eta*(y_xi*z_zeta-y_zeta*z_xi)   &
                            + x_zeta*(y_xi*z_eta-y_eta*z_xi)

                    det = 1.d0 / det

                    xi_x(i,j,k) = det*(y_eta*z_zeta - y_zeta*z_eta)
                    xi_y(i,j,k) = det*(x_zeta*z_eta - x_eta*z_zeta)
                    xi_z(i,j,k) = det*(x_eta*y_zeta - x_zeta*y_eta)

                    eta_x(i,j,k) = det*(y_zeta*z_xi - y_xi*z_zeta)
                    eta_y(i,j,k) = det*(x_xi*z_zeta - x_zeta*z_xi)
                    eta_z(i,j,k) = det*(x_zeta*y_xi - x_xi*y_zeta)

                    zeta_x(i,j,k) = det*(y_xi*z_eta - y_eta*z_xi)
                    zeta_y(i,j,k) = det*(x_eta*z_xi - x_xi*z_eta)
                    zeta_z(i,j,k) = det*(x_xi*y_eta - x_eta*y_xi)

                    !! Partial derivatives in global space
                    dfx_dx = fx_xi*xi_x(i,j,k) + fx_eta*eta_x(i,j,k) + fx_zeta*zeta_x(i,j,k)
                    dfy_dy = fy_xi*xi_y(i,j,k) + fy_eta*eta_y(i,j,k) + fy_zeta*zeta_y(i,j,k)
                    dfz_dz = fz_xi*xi_z(i,j,k) + fz_eta*eta_z(i,j,k) + fz_zeta*zeta_z(i,j,k)

                    !! Divergence
                    div_f(i,j,k) = dfx_dx + dfy_dy + dfz_dz

                end do
            end do
        end do

    end function divergence


    function gradient(f, x, y, z, imax, jmax, kmax) result(df)
        !! Finds the gradient of a scalar valued function f.
        implicit none
        integer, intent(in) :: imax, jmax, kmax
        real(8), dimension(imax,jmax,kmax), intent(in) :: f, x, y, z
        real(8), dimension(imax,jmax,kmax,3) :: df

        real(8) :: f_xi, f_eta, f_zeta
        real(8) :: x_xi, x_eta, x_zeta
        real(8) :: y_xi, y_eta, y_zeta
        real(8) :: z_xi, z_eta, z_zeta
        real(8), dimension(imax,jmax,kmax) :: xi_x, xi_y, xi_z
        real(8), dimension(imax,jmax,kmax) :: eta_x, eta_y, eta_z
        real(8), dimension(imax,jmax,kmax) :: zeta_x, zeta_y, zeta_z
        
        real(8) :: det
        integer :: i, j, k

        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax

                    !! Finite-Difference in parametric space.
                    if(i == 1) then
                        f_xi = f(2, j, k) - f(1, j, k)

                        x_xi = x(2, j, k) - x(1, j, k)
                        y_xi = y(2, j, k) - y(1, j, k)
                        z_xi = z(2, j, k) - z(1, j, k)
                    else if(i == imax) then
                        f_xi = f(imax, j, k) - f(imax-1, j, k)

                        x_xi = x(imax, j, k)-x(imax-1, j, k)
                        y_xi = y(imax, j, k)-y(imax-1, j, k)
                        z_xi = z(imax, j, k)-z(imax-1, j, k)
                    else
                        f_xi = 0.5d0*(f(i+1, j, k) - f(i-1, j, k))

                        x_xi = 0.5d0*(x(i+1, j, k) - x(i-1, j, k))
                        y_xi = 0.5d0*(y(i+1, j, k) - y(i-1, j, k))
                        z_xi = 0.5d0*(z(i+1, j, k) - z(i-1, j, k))
                    end if

                    if(j == 1) then
                        f_eta = f(i, 2, k) - f(i, 1, k)

                        x_eta = x(i, 2, k) - x(i, 1, k)
                        y_eta = y(i, 2, k) - y(i, 1, k)
                        z_eta = z(i, 2, k) - z(i, 1, k)
                    else if(j == jmax) then
                        f_eta = f(i, jmax, k) - f(i, jmax-1, k)

                        x_eta = x(i, jmax, k) - x(i, jmax-1, k)
                        y_eta = y(i, jmax, k) - y(i, jmax-1, k)
                        z_eta = z(i, jmax, k) - z(i, jmax-1, k)
                    else
                        f_eta = 0.5d0*(f(i, j+1, k) - f(i, j-1, k))

                        x_eta = 0.5d0*(x(i, j+1, k) - x(i, j-1, k))
                        y_eta = 0.5d0*(y(i, j+1, k) - y(i, j-1, k))
                        z_eta = 0.5d0*(z(i, j+1, k) - z(i, j-1, k))
                    end if

                    if(k == 1) then
                        f_zeta = f(i, j, 2) - f(i, j, 1)
                        x_zeta = x(i, j, 2) - x(i, j, 1)
                        y_zeta = y(i, j, 2) - y(i, j, 1)
                        z_zeta = z(i, j, 2) - z(i, j, 1)
                    else if(k == kmax) then
                        f_zeta = f(i, j, kmax) - f(i, j, kmax-1)

                        x_zeta = x(i, j, kmax) - x(i, j, kmax-1)
                        y_zeta = y(i, j, kmax) - y(i, j, kmax-1)
                        z_zeta = z(i, j, kmax) - z(i, j, kmax-1)
                    else
                        f_zeta = 0.5d0*(f(i, j, k+1) - f(i, j, k-1))

                        x_zeta = 0.5d0*(x(i, j, k+1) - x(i, j, k-1))
                        y_zeta = 0.5d0*(y(i, j, k+1) - y(i, j, k-1))
                        z_zeta = 0.5d0*(z(i, j, k+1) - z(i, j, k-1))
                    end if

                    !! Jacobian Transformation from parametric to global space.
                    det =   x_xi*(y_eta*z_zeta-y_zeta*z_eta)  &
                            - x_eta*(y_xi*z_zeta-y_zeta*z_xi)   &
                            + x_zeta*(y_xi*z_eta-y_eta*z_xi)

                    det = 1.d0 / det

                    xi_x(i,j,k) = det*(y_eta*z_zeta - y_zeta*z_eta)
                    xi_y(i,j,k) = det*(x_zeta*z_eta - x_eta*z_zeta)
                    xi_z(i,j,k) = det*(x_eta*y_zeta - x_zeta*y_eta)

                    eta_x(i,j,k) = det*(y_zeta*z_xi - y_xi*z_zeta)
                    eta_y(i,j,k) = det*(x_xi*z_zeta - x_zeta*z_xi)
                    eta_z(i,j,k) = det*(x_zeta*y_xi - x_xi*y_zeta)

                    zeta_x(i,j,k) = det*(y_xi*z_eta - y_eta*z_xi)
                    zeta_y(i,j,k) = det*(x_eta*z_xi - x_xi*z_eta)
                    zeta_z(i,j,k) = det*(x_xi*y_eta - x_eta*y_xi)

                    !! Partial derivatives in global space: (1) df_dx, (2) df_dy, (3) df_dz
                    df(i,j,k,1) = f_xi*xi_x(i,j,k) + f_eta*eta_x(i,j,k) + f_zeta*zeta_x(i,j,k)
                    df(i,j,k,2) = f_xi*xi_y(i,j,k) + f_eta*eta_y(i,j,k) + f_zeta*zeta_y(i,j,k)
                    df(i,j,k,3) = f_xi*xi_z(i,j,k) + f_eta*eta_z(i,j,k) + f_zeta*zeta_z(i,j,k)

                end do
            end do
        end do

    end function gradient

    
    function cross_product_3d(vec_1, vec_2) result(cross)
        !!! Cross product for 3 dimensional vectors.
        implicit none
        real(8), dimension(3), intent(in) :: vec_1, vec_2
        real(8), dimension(3) :: cross

        cross(1) = vec_1(2)*vec_2(3) - vec_1(3)*vec_2(2)
        cross(2) = vec_1(3)*vec_2(1) - vec_1(1)*vec_2(3)        
        cross(3) = vec_1(1)*vec_2(2) - vec_1(2)*vec_2(1)
    end function cross_product_3d

    
    function are_roots_real_cubic(a, b, c) result(are_real)
        !!! Determines if the roots of a cubic polynomical equation 
        !!! in the form of x^3 + a*x^2 + b*x + c = 0 are real.
        implicit none
        real(8), intent(in) :: a, b, c
        logical :: are_real
        real(8) :: q, r

        q = a*a - 3.d0*b / 9.d0
        r = (2.d0*a*a*a - 9.d0*a*b + 27.d0*c) / 54.d0

        are_real = (r*r < q*q*q)
    end function are_roots_real_cubic


    function roots_cubic_real(a, b, c) result(root)
        !!! Finds the real roots of a cubic polynomical equation 
        !!! in the form of x^3 + a*x^2 + b*x + c = 0.
        implicit none
        real(8), parameter :: pi = 4.d0*atan(1.d0)
        real(8), intent(in) :: a, b, c
        real(8), dimension(3) :: root
        real(8) :: q, r, theta

        q = a*a - 3.d0*b/9.d0
        r = (2.d0*a*a*a - 9.d0*a*b + 27.d0*c) / 54.d0
        theta = acos(r / sqrt(q*q*q))

        root(1) = -2.d0*sqrt(q)*cos(theta / 3.d0) - a / 3.d0
        root(2) = -2.d0*sqrt(q)*cos((theta + 2.d0*pi) / 3.d0) - a / 3.d0
        root(3) = -2.d0*sqrt(q)*cos((theta - 2.d0*pi) / 3.d0) - a / 3.d0
    end function roots_cubic_real


    function roots_cubic_imag(a, b, c) result(root)
        !!! Finds the imaginary roots of a cubic polynomical equation 
        !!! in the form of x^3 + a*x^2 + b*x + c = 0.
        implicit none
        real(8), parameter :: pi = 4.d0*atan(1.d0)
        complex(8), parameter :: i = (0.d0, 1.d0)   ! imaginary number i
        real(8), intent(in) :: a, b, c
        complex(8), dimension(3) :: root
        real(8) :: q, r, aa, bb

        q = a*a - 3.d0*b/9.d0
        r = (2.d0*a*a*a - 9.d0*a*b + 27.d0*c) / 54.d0

        aa = -sign(1.d0, r) * (abs(r) + sqrt(r*r - q*q*q))**(1.d0 / 3.d0)
        
        if (aa == 0.d0) then
            bb = 0.d0
        else
            bb = q / aa
        end if

        root(1) = aa + bb - a / 3.d0
        root(2) = cmplx(-0.5d0*(aa + bb) - (a / 3.d0), 0.5d0 * sqrt(3.d0) * (aa - bb), kind=8)
        root(3) = conjg(root(2))
    end function roots_cubic_imag


    function are_eig_real_3x3(mat) result(are_real)
        !!! Determines if the eigenvalues of a 3x3 matrix are real.
        real(8), dimension(3,3), intent(in) :: mat
        logical :: are_real
        real(8), dimension(3,3) :: t
        real(8) :: a, b, c
        
        a = -(mat(1,1) + mat(2,2) + mat(3,3))

        t = matmul(mat, mat)

        b = -0.5d0 * ( t(1,1) + t(2,2) + t(3,3) - (mat(1,1) + mat(2,2) + mat(3,3))**2 )

        c = -(mat(1,1) * (mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2))      &
            - mat(1,2) * (mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1))      &
            + mat(1,3) * (mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1)))

        are_real = are_roots_real_cubic(a, b, c)
    end function are_eig_real_3x3


    function eig_vals_real_3x3(mat) result(eig)
        !!! Finds the real eigenvalues of a 3x3 matrix
        implicit none
        real(8), dimension(3,3) :: mat
        complex(8), dimension(3) :: eig
        real(8), dimension(3,3) :: t
        real(8) :: a, b, c
        
        a = -(mat(1,1) + mat(2,2) + mat(3,3))

        t = matmul(mat, mat)

        b = -0.5d0 * ( t(1,1) + t(2,2) + t(3,3) - (mat(1,1) + mat(2,2) + mat(3,3))**2 )

        c = -(mat(1,1) * (mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2))      &
            - mat(1,2) * (mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1))      &
            + mat(1,3) * (mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1)))

        eig = roots_cubic_real(a, b, c)
    end function eig_vals_real_3x3


    function eig_vals_imag_3x3(mat) result(eig)
        !!! Finds the imaginary eigenvalues of 3x3 matrix.
        !!! eig(1) is real. 
        !!! eig(2) and eig(3) are imaginary.
        implicit none
        real(8), dimension(3,3) :: mat
        complex(8), dimension(3) :: eig
        real(8), dimension(3,3) :: t
        real(8) :: a, b, c
        
        a = -(mat(1,1) + mat(2,2) + mat(3,3))

        t = matmul(mat, mat)

        b = -0.5d0 * ( t(1,1) + t(2,2) + t(3,3) - (mat(1,1) + mat(2,2) + mat(3,3))**2 )

        c = -(mat(1,1) * (mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2))      &
            - mat(1,2) * (mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1))      &
            + mat(1,3) * (mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1)))
        
        eig = roots_cubic_imag(a, b, c)
    end function eig_vals_imag_3x3

    function negative_def(mat) result(neg_def)
        !!! Determines if a matrix is negative definite or not.
        implicit none
        real(8), dimension(3,3), intent(in) :: mat
        logical :: neg_def
        real(8) :: det
        logical :: check1, check2, check3

        check1 = mat(1,1) < 0.d0
        
        check2 = (mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)) > 0.d0

        det = mat(1,1)*(mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2))   &
            - mat(2,1)*(mat(1,2)*mat(3,3) - mat(1,3)*mat(3,2))   &
            + mat(3,1)*(mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2))

        check3 = det < 0.d0

        neg_def = (check1 .and. check2 .and. check3)
    end function negative_def

end module liutex_mod
