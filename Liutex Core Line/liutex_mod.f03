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


    function cross_product_3d(vec_1, vec_2) result(cross)
        !!! Cross product for 3 dimensional vectors.
        implicit none
        real(8), dimension(3), intent(in) :: vec_1, vec_2
        real(8), dimension(3) :: cross

        cross(1) = vec_1(2)*vec_2(3) - vec_1(3)*vec_2(2)
        cross(2) = vec_1(3)*vec_2(1) - vec_1(1)*vec_2(3)        
        cross(3) = vec_1(1)*vec_2(2) - vec_1(2)*vec_2(1)
    end function cross_product_3d


end module liutex_mod
