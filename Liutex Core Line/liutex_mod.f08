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


    !! Functions

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


    function finite_diff(f, x, y, z, imax, jmax, kmax) result(df)
        !! Finds the gradient of a scalar value f.
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

    end function finite_diff

    
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
