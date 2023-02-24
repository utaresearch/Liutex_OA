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
    function roots_cubic(a, b, c) result(root)
        !!! Finds the roots of a cubic polynomical equation 
        !!! in the form of x^3 + a*x^2 + b*x + c = 0.
        implicit none
        real(8), intent(in) :: a, b, c
        real(8), dimension(3) :: root

    end function roots_cubic

    function eig_vals(mat) result(eig)
        implicit none
        real(8), dimension(3,3) :: mat
        real(8), dimension(3) :: eig

    end function eig_vals

    function cross_product_3d(vec_1, vec_2) result(cross)
        !!! Cross product for 3 dimensional vectors.
        implicit none
        real(8), dimension(3), intent(in) :: vec_1, vec_2
        real(8), dimension(3) :: cross

        cross(1) = vec_1(2)*vec_2(3) - vec_1(3)*vec_2(2)
        cross(2) = vec_1(3)*vec_2(1) - vec_1(1)*vec_2(3)        
        cross(3) = vec_1(1)*vec_2(2) - vec_1(2)*vec_2(1)
    end function cross_product_3d

    
    function 
end module liutex_mod
