module liutex_UTA_mod
    !!-------------------------------------------------------------------------
    !! Fortran module that contains Liutex subroutine.
    !!
    !! By: Oscar Alvarez
    !! email: oscar.alvarez@uta.edu
    !! University of Texas at Arlington (UTA)
    !! Department of Mathematics
    !! Center for Numerical Simulation and Modeling (CNSM)
    !!-------------------------------------------------------------------------
    implicit none

    contains

    !! Subroutines
    subroutine liutex(a, liutex_vec, liutex_mag)
        !!! Calculates Liutex for a single point given its velocity gradient tensor.
        !!!
        !!! velocity gradient tensor := 
        !!! [ du/dx  du/dy  du/dz,
        !!!   dv/dx  dv/dy  dv/dz,
        !!!   dw/dx  dw/dy  dw/dz ].
        !!!

        implicit none

        real(8), dimension(3), parameter :: z0 = (/0.0, 0.0, 1.0/)

        real(8), dimension(3,3), intent(in) :: a            ! Velocity gradient tensor
        real(8), dimension(3), intent(out) :: liutex_vec    ! Liutex direction vector 
        real(8), intent(out) :: liutex_mag                  ! Liutex magnitude
      
        real(8), dimension(3,3) :: tt, qqq, vg
        real(8) :: aa, bb, cc, delta, rr, aaaa, bbbb, qq, delta1, delta2, delta3, temp
        real(8) :: alpha, beta
        real(8) :: eig3r
        complex(8) :: eig1c, eig2c 

        ! Cubic Formula
        ! Reference: Numerical Recipes in FORTRAN 77, Second Edition
        ! 5.6 Quadratic and Cubic Equations
        ! Page 179
        !---------------------------------------------------------------------

        ! cubic equation
        ! x**3 + aa * x**2 + bb * x + cc = 0

        ! coefficients of characteristic equation for the velocity gradient tensor.

        aa = -(a(1,1) + a(2,2) + a(3,3))

        tt = matmul(a,a)

        bb = -0.5d0*(tt(1,1) + tt(2,2) + tt(3,3) - (a(1,1) + a(2,2) + a(3,3))**2)

        cc = -( a(1,1) * (a(2,2)*a(3,3)-a(2,3)*a(3,2))                &
                - a(1,2) * (a(2,1)*a(3,3)-a(2,3)*a(3,1))              &
                + a(1,3) * (a(2,1)*a(3,2)-a(2,2)*a(3,1)) )

        ! discriminant of characteristic equation
        delta = 18.d0*aa*bb*cc - 4.d0*aa**3 * cc + aa**2 * bb**2 - 4.d0*bb**3 - 27.d0*cc**2

        ! alleviate round error
        delta = -delta / 108.d0
        liutex_vec = 0.d0
        liutex_mag = 0.d0

        if (delta > 0.d0) then ! one real root and two complex conjugate roots

            qq = (aa**2.d0 - 3.d0*bb) / 9.d0
            rr = (2.d0*aa**3 - 9.d0*aa*bb + 27.d0*cc) / 54.d0

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
                write(*,*) 'REAL EIG VALUE: ', eig3r
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
            
            call rotation(z0, liutex_vec, qqq)

            vg = matmul(transpose(qqq), a)
            vg = matmul(vg, qqq)

            alpha = 0.5d0 * sqrt((vg(2,2) - vg(1,1))**2 + (vg(2,1) + vg(1,2))**2)
            beta  = 0.5d0 * (vg(2,1) - vg(1,2))

            if(beta**2 > alpha**2) then

                if(beta > 0.0) then
                    liutex_mag = 2.d0 * (beta - alpha)
                    liutex_vec(1) = liutex_mag * liutex_vec(1)
                    liutex_vec(2) = liutex_mag * liutex_vec(2)
                    liutex_vec(3) = liutex_mag * liutex_vec(3)
                else
                    liutex_mag = 2.d0 * (beta + alpha)
                    liutex_vec(1) = liutex_mag * liutex_vec(1)
                    liutex_vec(2) = liutex_mag * liutex_vec(2)
                    liutex_vec(3) = liutex_mag * liutex_vec(3)
                end if

            else
                liutex_vec(1) = 0.d0
                liutex_vec(2) = 0.d0
                liutex_vec(3) = 0.d0
            end if

            liutex_mag = sqrt(liutex_vec(1)**2 + liutex_vec(2)**2 + liutex_vec(3)**2)
        else
            liutex_mag = 0.d0
            liutex_vec(1) = 0.d0
            liutex_vec(2) = 0.d0
            liutex_vec(3) = 0.d0
        end if

    end subroutine liutex


    subroutine rotation(u, v, r)
        !-------------------------------------------------------------------------------
        ! calculate rotation matrix r which rotates unit vector u to unit vector v
        !-------------------------------------------------------------------------------
        
        implicit none

        real(8), parameter :: eps = 1.0d-10
    
        real(8), dimension(3), intent(in) :: u, v
        real(8), dimension(3,3), intent(out) :: r
    
        real(8), dimension(3) :: a
        real(8) :: aa
        real(8) :: t
        real(8) :: alpha
        real(8) :: c, s

    
        ! a = u x v
        a(1) = u(2)*v(3) - u(3)*v(2)
        a(2) = u(3)*v(1) - u(1)*v(3)
        a(3) = u(1)*v(2) - u(2)*v(1)
    
        ! norm
        aa = sqrt(a(1)**2 + a(2)**2 + a(3)**2)
        
        if(aa < eps) then
            r = 0.d0
            
            r(1,1) = 1.d0
            r(2,2) = 1.d0
            r(3,3) = 1.d0
        else
            a = a/aa
            t = u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
            
            if(t > 1.d0) t = 1.d0
            if(t < -1.d0) t = -1.d0
            
            alpha = acos(t)
        
            c = cos(alpha)
            s = sin(alpha)
        
            r(1,1) = a(1)**2 * (1.d0 - c) + c
            r(1,2) = a(1) * a(2)*(1.d0 - c) - a(3) * s
            r(1,3) = a(1) * a(3) * (1.d0 - c) + a(2) * s
        
            r(2,1) = a(2) * a(1) * (1.d0 - c) + a(3)*s
            r(2,2) = a(2)**2 * (1.d0 - c) + c
            r(2,3) = a(2) * a(3) * (1.d0 - c) - a(1) * s
        
            r(3,1) = a(3) * a(1) * (1.d0 - c) - a(2) * s
            r(3,2) = a(3) * a(2) * (1.d0 - c) + a(1) * s
            r(3,3) = a(3)**2 * (1.d0 - c) + c
        end if
        
    end subroutine rotation

end module liutex_UTA_mod
