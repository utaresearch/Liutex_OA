module oa_stuff
contains

  subroutine finite_diff(du_dx, du_dy, du_dz, u, x, y, z, imax, jmax, kmax)
    !!! Calculates the partial derivatives dx, dy, dz using
    !!! up to 8th order finite difference.

    implicit none

    integer, intent(in) :: imax, jmax, kmax
    real, dimension(:,:,:), intent(in) :: x, y ,z, u
    real, dimension(:,:,:) :: du_dx, du_dy, du_dz
    
    real :: x_xi, x_eta, x_zeta
    real :: y_xi, y_eta, y_zeta
    real :: z_xi, z_eta, z_zeta
    real :: u_xi, u_eta, u_zeta

    real :: xi_x, eta_x, zeta_x
    real :: xi_y, eta_y, zeta_y
    real :: xi_z, eta_z, zeta_z
    real :: det

    integer :: i, j, k

    ! Finite-Difference for partial deriviates

    do k = 1, kmax
      do j = 1, jmax
        do i = 1, imax

          if(i == 1) then
            ! forward difference 6th order
            x_xi = ( -49.0*x(1,j,k) + 180.0*x(2,j,k) - 150.0*x(3,j,k) + 400.0*x(4,j,k)/3.0 &
                    - 75.0*x(5,j,k) + 24.0*x(6,j,k) - 10.0*x(7,j,k)/3.0 ) / 20.0
            y_xi = ( -49.0*y(1,j,k) + 180.0*y(2,j,k) - 150.0*y(3,j,k) + 400.0*y(4,j,k)/3.0 &
                    - 75.0*y(5,j,k) + 24.0*y(6,j,k) - 10.0*y(7,j,k)/3.0 ) / 20.0
            z_xi = ( -49.0*z(1,j,k) + 180.0*z(2,j,k) - 150.0*z(3,j,k) + 400.0*z(4,j,k)/3.0 &
                    - 75.0*z(5,j,k) + 24.0*z(6,j,k) - 10.0*z(7,j,k)/3.0 ) / 20.0

            u_xi = ( -49.0*u(1,j,k) + 180.0*u(2,j,k) - 150.0*u(3,j,k) + 400.0*u(4,j,k)/3.0 &
                    - 75.0*u(5,j,k) + 24.0*u(6,j,k) - 10.0*u(7,j,k)/3.0 ) / 20.0
          else if(i == imax) then
            ! backward difference 6th order
            x_xi = ( 49.0*x(imax,j,k) - 180.0*x(imax-1,j,k) + 150.0*x(imax-2,j,k) - 400.0*x(imax-3,j,k)/3.0 &
                    + 75.0*x(imax-4,j,k) - 24.0*x(imax-5,j,k) + 10.0*x(imax-6,j,k)/3.0 ) / 20.0
            y_xi = ( 49.0*y(imax,j,k) - 180.0*y(imax-1,j,k) + 150.0*y(imax-2,j,k) - 400.0*y(imax-3,j,k)/3.0 &
                    + 75.0*y(imax-4,j,k) - 24.0*y(imax-5,j,k) + 10.0*y(imax-6,j,k)/3.0 ) / 20.0
            z_xi = ( 49.0*z(imax,j,k) - 180.0*z(imax-1,j,k) + 150.0*z(imax-2,j,k) - 400.0*z(imax-3,j,k)/3.0 &
                    + 75.0*z(imax-4,j,k) - 24.0*z(imax-5,j,k) + 10.0*z(imax-6,j,k)/3.0 ) / 20.0

            u_xi = ( 49.0*u(imax,j,k) - 180.0*u(imax-1,j,k) + 150.0*u(imax-2,j,k) - 400.0*u(imax-3,j,k)/3.0 &
                    + 75.0*u(imax-4,j,k) - 24.0*u(imax-5,j,k) + 10.0*u(imax-6,j,k)/3.0 ) / 20.0
          else if(i == 2 .or. i == imax-1) then
            ! central difference 2nd order
            x_xi = 0.5*( x(i+1,j,k) - x(i-1,j,k) )
            y_xi = 0.5*( y(i+1,j,k) - y(i-1,j,k) )
            z_xi = 0.5*( z(i+1,j,k) - z(i-1,j,k) )

            u_xi = 0.5*( u(i+1,j,k) - u(i-1,j,k) )
          else if (i == 3 .or. i == imax-2) then
            ! central difference 4th order
            x_xi = ( -x(i+2,j,k) + 8.0*x(i+1,j,k) - 8.0*x(i-1,j,k) + x(i-2,j,k) ) / 12.0
            y_xi = ( -y(i+2,j,k) + 8.0*y(i+1,j,k) - 8.0*y(i-1,j,k) + y(i-2,j,k) ) / 12.0
            z_xi = ( -z(i+2,j,k) + 8.0*z(i+1,j,k) - 8.0*z(i-1,j,k) + z(i-2,j,k) ) / 12.0

            u_xi = ( -u(i+2,j,k) + 8.0*u(i+1,j,k) - 8.0*u(i-1,j,k) + u(i-2,j,k) ) / 12.0
          else if (i == 4 .or. i == imax-3) then
            ! central difference 6th order
            x_xi = ( -x(i-3,j,k) + 9.0*x(i-2,j,k) - 45.0*x(i-1,j,k) + 45.0*x(i+1,j,k) - 9.0*x(i+2,j,k) + x(i+3,j,k) ) / 60.0
            y_xi = ( -y(i-3,j,k) + 9.0*y(i-2,j,k) - 45.0*y(i-1,j,k) + 45.0*y(i+1,j,k) - 9.0*y(i+2,j,k) + y(i+3,j,k) ) / 60.0
            z_xi = ( -z(i-3,j,k) + 9.0*z(i-2,j,k) - 45.0*z(i-1,j,k) + 45.0*z(i+1,j,k) - 9.0*z(i+2,j,k) + z(i+3,j,k) ) / 60.0

            u_xi = ( -u(i-3,j,k) + 9.0*u(i-2,j,k) - 45.0*u(i-1,j,k) + 45.0*u(i+1,j,k) - 9.0*u(i+2,j,k) + u(i+3,j,k) ) / 60.0
          else
            ! central difference 8th order
            x_xi = ( x(i-4,j,k) - 32.0*x(i-3,j,k)/3.0 + 56.0*x(i-2,j,k) - 224.0*x(i-1,j,k) + 224.0*x(i+1,j,k) &
                    - 56.0*x(i+2,j,k) + 32.0*x(i+3,j,k)/3.0 - x(i+4,j,k) ) / 280.0
            y_xi = ( y(i-4,j,k) - 32.0*y(i-3,j,k)/3.0 + 56.0*y(i-2,j,k) - 224.0*y(i-1,j,k) + 224.0*y(i+1,j,k) &
                    - 56.0*y(i+2,j,k) + 32.0*y(i+3,j,k)/3.0 - y(i+4,j,k) ) / 280.0
            z_xi = ( z(i-4,j,k) - 32.0*z(i-3,j,k)/3.0 + 56.0*z(i-2,j,k) - 224.0*z(i-1,j,k) + 224.0*z(i+1,j,k) &
                    - 56.0*z(i+2,j,k) + 32.0*z(i+3,j,k)/3.0 - z(i+4,j,k) ) / 280.0

            u_xi = ( u(i-4,j,k) - 32.0*u(i-3,j,k)/3.0 + 56.0*u(i-2,j,k) - 224.0*u(i-1,j,k) + 224.0*u(i+1,j,k) &
                    - 56.0*u(i+2,j,k) + 32.0*u(i+3,j,k)/3.0 - u(i+4,j,k) ) / 280.0
          end if


          if(j == 1) then
            ! forward difference 6th order
            x_eta = ( -49.0*x(i,1,k) + 180.0*x(i,2,k) - 150.0*x(i,3,k) + 400.0*x(i,4,k)/3.0 &
                    - 75.0*x(i,5,k) + 24.0*x(i,6,k) - 10.0*x(i,7,k)/3.0 ) / 20.0
            y_eta = ( -49.0*y(i,1,k) + 180.0*y(i,2,k) - 150.0*y(i,3,k) + 400.0*y(i,4,k)/3.0 &
                    - 75.0*y(i,5,k) + 24.0*y(i,6,k) - 10.0*y(i,7,k)/3.0 ) / 20.0
            z_eta = ( -49.0*z(i,1,k) + 180.0*z(i,2,k) - 150.0*z(i,3,k) + 400.0*z(i,4,k)/3.0 &
                    - 75.0*z(i,5,k) + 24.0*z(i,6,k) - 10.0*z(i,7,k)/3.0 ) / 20.0

            u_eta = ( -49.0*u(i,1,k) + 180.0*u(i,2,k) - 150.0*u(i,3,k) + 400.0*u(i,4,k)/3.0 &
                    - 75.0*u(i,5,k) + 24.0*u(i,6,k) - 10.0*u(i,7,k)/3.0 ) / 20.0
          else if(j == jmax) then
            ! backward difference 6th order
            x_eta = ( 49.0*x(i,jmax,k) - 180.0*x(i,jmax-1,k) + 150.0*x(i,jmax-2,k) - 400.0*x(i,jmax-3,k)/3.0 &
                    + 75.0*x(i,jmax-4,k) - 24.0*x(i,jmax-5,k) + 10.0*x(i,jmax-6,k)/3.0 ) / 20.0
            y_eta = ( 49.0*y(i,jmax,k) - 180.0*y(i,jmax-1,k) + 150.0*y(i,jmax-2,k) - 400.0*y(i,jmax-3,k)/3.0 &
                    + 75.0*y(i,jmax-4,k) - 24.0*y(i,jmax-5,k) + 10.0*y(i,jmax-6,k)/3.0 ) / 20.0
            z_eta = ( 49.0*z(i,jmax,k) - 180.0*z(i,jmax-1,k) + 150.0*z(i,jmax-2,k) - 400.0*z(i,jmax-3,k)/3.0 &
                    + 75.0*z(i,jmax-4,k) - 24.0*z(i,jmax-5,k) + 10.0*z(i,jmax-6,k)/3.0 ) / 20.0

            u_eta = ( 49.0*u(i,jmax,k) - 180.0*u(i,jmax-1,k) + 150.0*u(i,jmax-2,k) - 400.0*u(i,jmax-3,k)/3.0 &
                    + 75.0*u(i,jmax-4,k) - 24.0*u(i,jmax-5,k) + 10.0*u(i,jmax-6,k)/3.0 ) / 20.0
          else if(j == 2 .or. j == jmax-1) then
            ! central difference 2nd order
            x_eta = 0.5*( x(i,j+1,k) - x(i,j-1,k) )
            y_eta = 0.5*( y(i,j+1,k) - y(i,j-1,k) )
            z_eta = 0.5*( z(i,j+1,k) - z(i,j-1,k) )

            u_eta = 0.5*( u(i,j+1,k) - u(i,j-1,k) )
          else if (j == 3 .or. j == jmax-2) then
            ! central difference 4th order
            x_eta = ( -x(i,j+2,k) + 8.0*x(i,j+1,k) - 8.0*x(i,j-1,k) + x(i,j-2,k) )/12.0
            y_eta = ( -y(i,j+2,k) + 8.0*y(i,j+1,k) - 8.0*y(i,j-1,k) + y(i,j-2,k) )/12.0
            z_eta = ( -z(i,j+2,k) + 8.0*z(i,j+1,k) - 8.0*z(i,j-1,k) + z(i,j-2,k) )/12.0

            u_eta = ( -u(i,j+2,k) + 8.0*u(i,j+1,k) - 8.0*u(i,j-1,k) + u(i,j-2,k) )/12.0
          else if (j == 4 .or. j == jmax-3) then
            ! central difference 6th order
            x_eta = ( -x(i,j-3,k) + 9.0*x(i,j-2,k) - 45.0*x(i,j-1,k) + 45.0*x(i,j+1,k) - 9.0*x(i,j+2,k) + x(i,j+3,k) ) / 60.0
            y_eta = ( -y(i,j-3,k) + 9.0*y(i,j-2,k) - 45.0*y(i,j-1,k) + 45.0*y(i,j+1,k) - 9.0*y(i,j+2,k) + y(i,j+3,k) ) / 60.0
            z_eta = ( -z(i,j-3,k) + 9.0*z(i,j-2,k) - 45.0*z(i,j-1,k) + 45.0*z(i,j+1,k) - 9.0*z(i,j+2,k) + z(i,j+3,k) ) / 60.0

            u_eta = ( -u(i,j-3,k) + 9.0*u(i,j-2,k) - 45.0*u(i,j-1,k) + 45.0*u(i,j+1,k) - 9.0*u(i,j+2,k) + u(i,j+3,k) ) / 60.0
          else
            ! central difference 8th order
            x_eta = ( x(i,j-4,k) - 32.0*x(i,j-3,k)/3.0 + 56.0*x(i,j-2,k) - 224.0*x(i,j-1,k) + 224.0*x(i,j+1,k) &
                    - 56.0*x(i,j+2,k) + 32.0*x(i,j+3,k)/3.0 - x(i,j+4,k) ) / 280.0
            y_eta = ( y(i,j-4,k) - 32.0*y(i,j-3,k)/3.0 + 56.0*y(i,j-2,k) - 224.0*y(i,j-1,k) + 224.0*y(i,j+1,k) &
                    - 56.0*y(i,j+2,k) + 32.0*y(i,j+3,k)/3.0 - y(i,j+4,k) ) / 280.0
            z_eta = ( z(i,j-4,k) - 32.0*z(i,j-3,k)/3.0 + 56.0*z(i,j-2,k) - 224.0*z(i,j-1,k) + 224.0*z(i,j+1,k) &
                    - 56.0*z(i,j+2,k) + 32.0*z(i,j+3,k)/3.0 - z(i,j+4,k) ) / 280.0

            u_eta = ( u(i,j-4,k) - 32.0*u(i,j-3,k)/3.0 + 56.0*u(i,j-2,k) - 224.0*u(i,j-1,k) + 224.0*u(i,j+1,k) &
                    - 56.0*u(i,j+2,k) + 32.0*u(i,j+3,k)/3.0 - u(i,j+4,k) ) / 280.0
          end if


          if(k == 1) then
            ! forward difference 6th order
            x_zeta = ( -49.0*x(i,j,1) + 180.0*x(i,j,2) - 150.0*x(i,j,3) + 400.0*x(i,j,4)/3.0 &
                    - 75.0*x(i,j,5) + 24.0*x(i,j,6) - 10.0*x(i,j,7)/3.0 ) / 20.0
            y_zeta = ( -49.0*y(i,j,1) + 180.0*y(i,j,2) - 150.0*y(i,j,3) + 400.0*y(i,j,4)/3.0 &
                    - 75.0*y(i,j,5) + 24.0*y(i,j,6) - 10.0*y(i,j,7)/3.0 ) / 20.0
            z_zeta = ( -49.0*z(i,j,1) + 180.0*z(i,j,2) - 150.0*z(i,j,3) + 400.0*z(i,j,4)/3.0 &
                    - 75.0*z(i,j,5) + 24.0*z(i,j,6) - 10.0*z(i,j,7)/3.0 ) / 20.0

            u_zeta = ( -49.0*u(i,j,1) + 180.0*u(i,j,2) - 150.0*u(i,j,3) + 400.0*u(i,j,4)/3.0 &
                    - 75.0*u(i,j,5) + 24.0*u(i,j,6) - 10.0*u(i,j,7)/3.0 ) / 20.0
          else if(k == kmax) then
            ! backward difference 6th order
            x_zeta = ( 49.0*x(i,j,kmax) - 180.0*x(i,j,kmax-1) + 150.0*x(i,j,kmax-2) - 400.0*x(i,j,kmax-3)/3.0 &
                    + 75.0*x(i,j,kmax-4) - 24.0*x(i,j,kmax-5) + 10.0*x(i,j,kmax-6)/3.0 ) / 20.0
            y_zeta = ( 49.0*y(i,j,kmax) - 180.0*y(i,j,kmax-1) + 150.0*y(i,j,kmax-2) - 400.0*y(i,j,kmax-3)/3.0 &
                    + 75.0*y(i,j,kmax-4) - 24.0*y(i,j,kmax-5) + 10.0*y(i,j,kmax-6)/3.0 ) / 20.0
            z_zeta = ( 49.0*z(i,j,kmax) - 180.0*z(i,j,kmax-1) + 150.0*z(i,j,kmax-2) - 400.0*z(i,j,kmax-3)/3.0 &
                    + 75.0*z(i,j,kmax-4) - 24.0*z(i,j,kmax-5) + 10.0*z(i,j,kmax-6)/3.0 ) / 20.0

            u_zeta = ( 49.0*u(i,j,kmax) - 180.0*u(i,j,kmax-1) + 150.0*u(i,j,kmax-2) - 400.0*u(i,j,kmax-3)/3.0 &
                    + 75.0*u(i,j,kmax-4) - 24.0*u(i,j,kmax-5) + 10.0*u(i,j,kmax-6)/3.0 ) / 20.0
          else if(k == 2 .or. k == kmax-1) then
            ! central difference 2nd order
            x_zeta = 0.5*( x(i,j,k+1) - x(i,j,k-1) )
            y_zeta = 0.5*( y(i,j,k+1) - y(i,j,k-1) )
            z_zeta = 0.5*( z(i,j,k+1) - z(i,j,k-1) )

            u_zeta = 0.5*( u(i,j,k+1) - u(i,j,k-1) )
          else if (k == 3 .or. k == kmax-2) then
            ! central difference
            x_zeta = ( -x(i,j,k+2) + 8.0*x(i,j,k+1) - 8.0*x(i,j,k-1) + x(i,j,k-2) )/12.0
            y_zeta = ( -y(i,j,k+2) + 8.0*y(i,j,k+1) - 8.0*y(i,j,k-1) + y(i,j,k-2) )/12.0
            z_zeta = ( -z(i,j,k+2) + 8.0*z(i,j,k+1) - 8.0*z(i,j,k-1) + z(i,j,k-2) )/12.0

            u_zeta = ( -u(i,j,k+2) + 8.0*u(i,j,k+1) - 8.0*u(i,j,k-1) + u(i,j,k-2) )/12.0
          else if (k == 4 .or. k == kmax-3) then
            ! central difference 6th order
            x_zeta = ( -x(i,j,k-3) + 9.0*x(i,j,k-2) - 45.0*x(i,j,k-1) + 45.0*x(i,j,k+1) - 9.0*x(i,j,k+2) + x(i,j,k+3) ) / 60.0
            y_zeta = ( -y(i,j,k-3) + 9.0*y(i,j,k-2) - 45.0*y(i,j,k-1) + 45.0*y(i,j,k+1) - 9.0*y(i,j,k+2) + y(i,j,k+3) ) / 60.0
            z_zeta = ( -z(i,j,k-3) + 9.0*z(i,j,k-2) - 45.0*z(i,j,k-1) + 45.0*z(i,j,k+1) - 9.0*z(i,j,k+2) + z(i,j,k+3) ) / 60.0

            u_zeta = ( -u(i,j,k-3) + 9.0*u(i,j,k-2) - 45.0*u(i,j,k-1) + 45.0*u(i,j,k+1) - 9.0*u(i,j,k+2) + u(i,j,k+3) ) / 60.0
          else
            ! central difference 8th order
            x_zeta = ( x(i,j,k-4) - 32.0*x(i,j,k-3)/3.0 + 56.0*x(i,j,k-2) - 224.0*x(i,j,k-1) + 224.0*x(i,j,k+1) &
                    - 56.0*x(i,j,k+2) + 32.0*x(i,j,k+3)/3.0 - x(i,j,k+4) ) / 280.0
            y_zeta = ( y(i,j,k-4) - 32.0*y(i,j,k-3)/3.0 + 56.0*y(i,j,k-2) - 224.0*y(i,j,k-1) + 224.0*y(i,j,k+1) &
                    - 56.0*y(i,j,k+2) + 32.0*y(i,j,k+3)/3.0 - y(i,j,k+4) ) / 280.0
            z_zeta = ( z(i,j,k-4) - 32.0*z(i,j,k-3)/3.0 + 56.0*z(i,j,k-2) - 224.0*z(i,j,k-1) + 224.0*z(i,j,k+1) &
                    - 56.0*z(i,j,k+2) + 32.0*z(i,j,k+3)/3.0 - z(i,j,k+4) ) / 280.0

            u_zeta = ( u(i,j,k-4) - 32.0*u(i,j,k-3)/3.0 + 56.0*u(i,j,k-2) - 224.0*u(i,j,k-1) + 224.0*u(i,j,k+1) &
                    - 56.0*u(i,j,k+2) + 32.0*u(i,j,k+3)/3.0 - u(i,j,k+4) ) / 280.0
          end if

          ! determinant of Jacobian
          det =   x_xi*(y_eta*z_zeta-y_zeta*z_eta)  &
                  - x_eta*(y_xi*z_zeta-y_zeta*z_xi)   &
                  + x_zeta*(y_xi*z_eta-y_eta*z_xi)

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

          ! Liutex Magnitude Gradient (first partial derivatives)

          du_dx(i,j,k) = u_xi*xi_x + u_eta*eta_x + u_zeta*zeta_x
          du_dy(i,j,k) = u_xi*xi_y + u_eta*eta_y + u_zeta*zeta_y
          du_dz(i,j,k) = u_xi*xi_z + u_eta*eta_z + u_zeta*zeta_z

        end do
      end do
    end do

  end subroutine

end module oa_stuff


program DNSUTA_extract
!-------------------------------------------------------------------------------
! Author: Oscar Alvarez
! Email: oscar.alvarez@mavs.uta.edu
! Email: gaoyisheng@me.com, cliu@uta.edu
! Department of Mathematics, University of Texas at Arlington,
! Arlington, Texas, USA
! Fall 2022
!-------------------------------------------------------------------------------
  use oa_stuff

  implicit none

  integer :: imax, jmax, kmax
  integer :: imax_local, jmax_local, kmax_local

  integer :: istart, iend
  integer :: jstart, jend
  integer :: kstart, kend

  real, dimension(:,:,:), allocatable :: x_global, x
  real, dimension(:,:,:), allocatable :: y_global, y
  real, dimension(:,:,:), allocatable :: z_global, z

  real, dimension(:,:,:), allocatable :: u
  real, dimension(:,:,:), allocatable :: v
  real, dimension(:,:,:), allocatable :: w

  real, dimension(:,:,:,:), allocatable :: f

  real, dimension(:,:,:), allocatable :: dudx, dudy, dudz
  real, dimension(:,:,:), allocatable :: dvdx, dvdy, dvdz
  real, dimension(:,:,:), allocatable :: dwdx, dwdy, dwdz
  real :: rdudx,rdudy,rdudz,rdvdx,rdvdy,rdvdz,rdwdx
  real :: rdwdy,rdwdz,rvort_x,rvort_y,rvort_z

  real, dimension(:,:,:), allocatable :: vorticity_x
  real, dimension(:,:,:), allocatable :: vorticity_y
  real, dimension(:,:,:), allocatable :: vorticity_z

  real, dimension(:,:,:), allocatable :: vorticity_mag

  real, dimension(:,:,:), allocatable :: omega
  real, dimension(:,:,:), allocatable :: omega_x, omega_y, omega_z, omgGradXliutex
  real, dimension(:,:,:), allocatable :: omega_xx, omega_xy, omega_xz
  real, dimension(:,:,:), allocatable :: omega_yx, omega_yy, omega_yz
  real, dimension(:,:,:), allocatable :: omega_zx, omega_zy, omega_zz

  real, dimension(:,:,:), allocatable :: q, q_x, q_y, q_z

  real, dimension(:,:,:), allocatable :: liutex_x
  real, dimension(:,:,:), allocatable :: liutex_y
  real, dimension(:,:,:), allocatable :: liutex_z
  real, dimension(:,:,:), allocatable :: isLocExtrOmg
  real, dimension(:,:,:), allocatable :: liutexEigr

  real, dimension(:,:,:), allocatable :: liutex_mag, lamb_ci
  real,dimension (:,:,:), allocatable :: liutex_mag_x, liutex_mag_y, liutex_mag_z

  real, dimension(3) :: vor, vrtmp
  real :: rorMag, ljmtmp, ljmtmp2

  character(100) :: gridfilename, funcfilename, outputfilename
  character(100) :: outputfilename1, outputfilename2, outputfilename3
  integer :: skip
  integer :: f_start
  integer :: f_end

  integer :: iii

  integer :: i, j, k, nx
  integer :: ii, jj, kk

  integer, parameter :: fin1 = 10
  integer, parameter :: fin2 = 11

  integer, parameter :: fout = 12

  logical :: lexist

  integer :: ios
  character(100) :: msg
  integer :: nvar

  real :: aaa
  real :: bbb

  real :: thljm, qljm,rljm,pi
  real :: omega_eps

  real :: vort_2

  real :: a(3,3), Aljm(3,3), Matljm(3,3)

  real :: t1, t2, t5, t6

  real :: aa, bb, cc
  real :: delta
  real :: tt(3, 3)

  complex :: eig1c, eig2c
  real :: eig3r

  real :: qq, rr, root11,root22,root33
  real :: aaaa, bbbb

  real :: vr(3)
  real :: temp
  real :: vg(3,3)

  real :: alpha, beta

  ! Modified Omega Liutex variables
  real, dimension(:,:,:), allocatable :: omega_l
  real, dimension(:,:,:), allocatable :: lambda_cr, lambda_r
  real :: maxbeta_alpha, omega_l_eps

  real, dimension(:,:,:), allocatable :: o_alpha, o_beta

  real, parameter :: z0(3) = (/0.0, 0.0, 1.0/)
  real :: rm
  real :: qqq(3,3)

  real :: delta1, delta2, delta3

  character(6) :: chars

  ! the default epsilon value for Omega method
  omega_eps = 1.0e-4
  pi = 4.0 * atan(1.0)

  ! inquire whether the input file exists
  inquire(file='input.txt', exist=lexist)
  if(.not. lexist) then
    write(*,*) 'ERROR: no file input.txt'
    stop
  end if

  ! read input file
  open(fin1, file='input.txt', form='formatted', action='read')
  read(fin1, *) f_start
  read(fin1, *) f_end
  read(fin1, *) skip
  read(fin1, *) istart, iend
  read(fin1, *) jstart, jend
  read(fin1, *) kstart, kend
  read(fin1, *) omega_eps
  read(fin1, *) outputfilename
  close(fin1)

  imax_local = iend - istart + 1
  jmax_local = jend - jstart + 1
  kmax_local = kend - kstart + 1

  omega_l_eps = omega_eps

  gridfilename = 'data/plate.xyz'

  ! check if the grid file exists
  inquire(file=trim(gridfilename), exist=lexist)
  if(.not. lexist) then
    write(*,*) 'ERROR: no grid file ', trim(gridfilename)
    stop
  end if

  write(*,*)
  write(*,*) 'reading ', trim(gridfilename)

  ! open grid file
  open(fin1, file=trim(gridfilename), form='unformatted', action='read',       &
       iostat=ios, iomsg=msg)
  if(ios /= 0) then
    write(*,*) 'ERROR: ', msg
    stop
  end if

  read(fin1, iostat=ios, iomsg=msg) imax, jmax, kmax
  if(ios /= 0) then
    write(*,*) 'ERROR: ', msg
    stop
  end if

  allocate(x_global(imax, jmax, kmax))
  allocate(y_global(imax, jmax, kmax))
  allocate(z_global(imax, jmax, kmax))

  write(*,*)
  write(*,*) 'imax = ', imax, ' jmax = ', jmax, ' kmax = ', kmax

  ! read coordinates
  read(fin1, iostat=ios, iomsg=msg)                                            &
                              (((x_global(i, j, k), i=1,imax), j=1,jmax), k=1,kmax),  &
                              (((y_global(i, j, k), i=1,imax), j=1,jmax), k=1,kmax),  &
                              (((z_global(i, j, k), i=1,imax), j=1,jmax), k=1,kmax)

  close(fin1)


  do iii = f_start, f_end, skip

    write(*,*) 'imax_local, jmax_local, kmax_local'
    write(*,*) imax_local, jmax_local, kmax_local

    call itoa6(iii, chars)

    funcfilename = 'data/plate_'//chars//'.fun'

    ! check if the data file exists
    inquire(file=trim(funcfilename), exist=lexist)
    if(.not. lexist) then
      write(*,*) 'ERROR: no data file ', funcfilename
      stop
    end if

    write(*,*)
    write(*,*) 'reading ', trim(funcfilename)

    ! open data file
    open(fin2, file=trim(funcfilename), form='unformatted', action='read',     &
         iostat=ios, iomsg=msg)
    if(ios /= 0) then
      write(*,*) 'ERROR: ', msg
      stop
    end if

    read(fin2, iostat=ios, iomsg=msg) imax, jmax, kmax, nvar
    if(ios /= 0) then
      write(*,*) 'ERROR: ', msg
      stop
    end if

    allocate(f(imax, jmax, kmax, nvar))

    ! read data
    !---------------------------------------------------------------------------
    ! f(1): density
    ! f(2): u velocity
    ! f(3): v velocity
    ! f(4): w velocity
    ! f(5): pressure
    !---------------------------------------------------------------------------

    read(fin2, iostat=ios, iomsg=msg) ((((f(i, j, k, nx), i=1,imax), j=1,jmax),&
                                      k=1,kmax), nx=1,nvar)
    if(ios /= 0) then
      write(*,*) 'ERROR: ', msg
      stop
    end if

    close(fin2)

    allocate(dudx(imax_local, jmax_local, kmax_local))
    allocate(dudy(imax_local, jmax_local, kmax_local))
    allocate(dudz(imax_local, jmax_local, kmax_local))
    allocate(dvdx(imax_local, jmax_local, kmax_local))
    allocate(dvdy(imax_local, jmax_local, kmax_local))
    allocate(dvdz(imax_local, jmax_local, kmax_local))
    allocate(dwdx(imax_local, jmax_local, kmax_local))
    allocate(dwdy(imax_local, jmax_local, kmax_local))
    allocate(dwdz(imax_local, jmax_local, kmax_local))

    allocate(vorticity_x(imax_local,jmax_local,kmax_local))
    allocate(vorticity_y(imax_local,jmax_local,kmax_local))
    allocate(vorticity_z(imax_local,jmax_local,kmax_local))

    allocate(vorticity_mag(imax_local,jmax_local,kmax_local))

    allocate(omega(imax_local, jmax_local, kmax_local))

    allocate(q(imax_local, jmax_local, kmax_local))

    allocate(q_x(imax_local, jmax_local, kmax_local))
    allocate(q_y(imax_local, jmax_local, kmax_local))
    allocate(q_z(imax_local, jmax_local, kmax_local))

    allocate(omgGradXliutex(imax_local, jmax_local, kmax_local))

    allocate(omega_x(imax_local, jmax_local, kmax_local))
    allocate(omega_y(imax_local, jmax_local, kmax_local))
    allocate(omega_z(imax_local, jmax_local, kmax_local))

    allocate(omega_xx(imax_local, jmax_local, kmax_local))
    allocate(omega_xy(imax_local, jmax_local, kmax_local))
    allocate(omega_xz(imax_local, jmax_local, kmax_local))

    allocate(omega_yx(imax_local, jmax_local, kmax_local))
    allocate(omega_yy(imax_local, jmax_local, kmax_local))
    allocate(omega_yz(imax_local, jmax_local, kmax_local))

    allocate(omega_zx(imax_local, jmax_local, kmax_local))
    allocate(omega_zy(imax_local, jmax_local, kmax_local))
    allocate(omega_zz(imax_local, jmax_local, kmax_local))

    allocate(isLocExtrOmg(imax_local, jmax_local, kmax_local))

    allocate(u(imax_local, jmax_local, kmax_local))
    allocate(v(imax_local, jmax_local, kmax_local))
    allocate(w(imax_local, jmax_local, kmax_local))

    allocate(liutex_x(imax_local, jmax_local, kmax_local))
    allocate(liutex_y(imax_local, jmax_local, kmax_local))
    allocate(liutex_z(imax_local, jmax_local, kmax_local))

    allocate(liutexEigr(imax_local, jmax_local, kmax_local))

    allocate(liutex_mag(imax_local, jmax_local, kmax_local))
    allocate(liutex_mag_x(imax_local, jmax_local, kmax_local))
    allocate(liutex_mag_y(imax_local, jmax_local, kmax_local))
    allocate(liutex_mag_z(imax_local, jmax_local, kmax_local))

    allocate(lamb_ci(imax_local, jmax_local, kmax_local))

    isLocExtrOmg = 0.0
    omgGradXliutex = 0.0
    liutexEigr = 0.0

    dudx = 0.0
    dudy = 0.0
    dudz = 0.0
    dvdx = 0.0
    dvdy = 0.0
    dvdz = 0.0
    dwdx = 0.0
    dwdy = 0.0
    dwdz = 0.0

    omega = 0.0
    q = 0.0

    !---------------------------------------------------------------------------
    ! calculate:
    ! (1) velocity gradient tensor
    ! (2) vorticity
    ! (3) Omega
    ! (4) Q
    ! For coordinate transformation, refer to CFL3D User's Manual Appendix F
    !---------------------------------------------------------------------------

    write(*,*)
    write(*,*) 'Calculating velocity gradient tensor'

    call cpu_time(t1)

    allocate(x(imax_local, jmax_local, kmax_local))
    allocate(y(imax_local, jmax_local, kmax_local))
    allocate(z(imax_local, jmax_local, kmax_local))

    do k = kstart, kend
      do j = jstart, jend
        do i = istart, iend

          ii = i - istart + 1
          jj = j - jstart + 1
          kk = k - kstart + 1

!          write(*,*) 'ii = ', ii, 'jj = ', jj, 'kk = ', kk

          x(ii, jj, kk) = x_global(i, j, k)
          y(ii, jj, kk) = y_global(i, j, k)
          z(ii, jj, kk) = z_global(i ,j, k)

!          write(*,*) 'xyz passed'

          u(ii, jj, kk) = f(i, j, k, 2)
          v(ii, jj, kk) = f(i, j, k, 3)
          w(ii, jj, kk) = f(i ,j, k, 4)

!          write(*,*) 'uvw passed'
        end do
      end do
    end do

    write(*,*) 'passed variable index shift'

    call finite_diff(dudx, dudy, dudz, u, x, y, z, imax_local, jmax_local, kmax_local)
    call finite_diff(dvdx, dvdy, dvdz, v, x, y, z, imax_local, jmax_local, kmax_local)
    call finite_diff(dwdx, dwdy, dwdz, w, x, y, z, imax_local, jmax_local, kmax_local)

    !--------

    do k = 1, kmax_local
      do j = 1, jmax_local
        do i = 1, imax_local

          rdudx = dudx(i,j,k)
          rdudy = dudy(i,j,k)
          rdudz = dudz(i,j,k)

          rdvdx = dvdx(i,j,k)
          rdvdy = dvdy(i,j,k)
          rdvdz = dvdz(i,j,k)

          rdwdx = dwdx(i,j,k)
          rdwdy = dwdy(i,j,k)
          rdwdz = dwdz(i,j,k)

          aaa = rdudx**2 + rdvdy**2 + rdwdz**2 +               &
                0.5*((rdudy+rdvdx)**2 +                        &
                     (rdudz+rdwdx)**2 +                        &
                     (rdvdz+rdwdy)**2)

          ! vorticity
          rvort_x = rdwdy - rdvdz
          rvort_y = rdudz - rdwdx
          rvort_z = rdvdx - rdudy

          vor(1) = rvort_x
          vor(2) = rvort_y
          vor(3) = rvort_z

          vort_2 = rvort_x**2 + rvort_y**2 + rvort_z**2

          bbb = 0.5 * vort_2

          omega(i,j,k) = bbb / (aaa + bbb + omega_eps)

          a(1,1) = rdudx
          a(1,2) = rdudy
          a(1,3) = rdudz
          a(2,1) = rdvdx
          a(2,2) = rdvdy
          a(2,3) = rdvdz
          a(3,1) = rdwdx
          a(3,2) = rdwdy
          a(3,3) = rdwdz

          call cal_liutex(a, vor, vrtmp, rorMag)

          liutex_mag(i,j,k) = rorMag

          q(i,j,k) = 0.5 * (bbb - aaa)

        end do
      end do
    end do

    !+++++++

    call cpu_time(t2)

    write(*,*) 'calculation time: ', t2 - t1

    ! calculate hessain matrix of the omega
    ! fist step to calculate the first order dirivative of the omega

    call finite_diff(omega_x, omega_y, omega_z, omega, x, y, z, imax_local, jmax_local, kmax_local)
    call finite_diff(liutex_mag_x, liutex_mag_y, liutex_mag_z, liutex_mag, x, y, z, imax_local, jmax_local, kmax_local)
    call finite_diff(q_x, q_y, q_z, q, x, y, z, imax_local, jmax_local, kmax_local)

    ! end of the first step of hessain matrix computation

    call finite_diff(omega_xx, omega_xy, omega_xz, omega_x, x, y, z, imax_local, jmax_local, kmax_local)
    call finite_diff(omega_yx, omega_yy, omega_yz, omega_y, x, y, z, imax_local, jmax_local, kmax_local)
    call finite_diff(omega_zx, omega_zy, omega_zz, omega_z, x, y, z, imax_local, jmax_local, kmax_local)

    !calculate the eigenvalue

    do k = 1, kmax_local
      do j = 1, jmax_local
        do i = 1, imax_local

          ! set velocity gradient tensor
          a(1,1) = omega_xx(i,j,k)
          a(1,2) = omega_xy(i,j,k)
          a(1,3) = omega_xz(i,j,k)
          a(2,1) = omega_yx(i,j,k)
          a(2,2) = omega_yy(i,j,k)
          a(2,3) = omega_yz(i,j,k)
          a(3,1) = omega_zx(i,j,k)
          a(3,2) = omega_zy(i,j,k)
          a(3,3) = omega_zz(i,j,k)

          Aljm = (a + transpose(a)) / 2

          Matljm = Aljm
          aa = -(Matljm(1,1)+Matljm(2,2) + Matljm(3,3))

          tt = matmul(Matljm,Matljm)

          bb = -0.5*(tt(1,1) + tt(2,2) + tt(3,3) - (Matljm(1,1) + Matljm(2,2) + Matljm(3,3))**2)

          cc = -(Matljm(1,1)*(Matljm(2,2)*Matljm(3,3)-Matljm(2,3)*Matljm(3,2))                            &
                 -Matljm(1,2)*(Matljm(2,1)*Matljm(3,3)-Matljm(2,3)*Matljm(3,1))                           &
                 +Matljm(1,3)*(Matljm(2,1)*Matljm(3,2)-Matljm(2,2)*Matljm(3,1)))

          qljm = (aa**2-3*bb)/9
          rljm = (2*aa**3-9*aa*bb+27*cc)/54
          thljm = acos(rljm/sqrt(qljm**3))
          root11 = -2*sqrt(qljm)*cos(thljm/3)-aa/3
          root22 = -2*sqrt(qljm)*cos(thljm/3+2.0/3*pi)-aa/3
          root33 = -2*sqrt(qljm)*cos(thljm/3-2.0/3*pi)-aa/3

          if (root11>root22) then
              qljm = root11
              root11 = root22
              root22 = qljm
          end if

          if (root22 > root33) then
              qljm = root22
              root22 = root33
              root33 = qljm
          end if

          if (root11 > root22) then
              qljm = root11
              root11 = root22
              root22 = qljm
          end if

          if (root22 < 0.0) then
            isLocExtrOmg(i,j,k) = 1.0
          end if

        end do
      end do
    end do

    deallocate(omega_x)
    deallocate(omega_y)
    deallocate(omega_z)

    deallocate(omega_xx)
    deallocate(omega_xy)
    deallocate(omega_xz)

    deallocate(omega_yx)
    deallocate(omega_yy)
    deallocate(omega_yz)

    deallocate(omega_zx)
    deallocate(omega_zy)
    deallocate(omega_zz)

    ! end of the hessain matrix computation

    deallocate(f)

    !---------------------------------------------------------------------------
    ! calculate liutex
    !---------------------------------------------------------------------------

    lamb_ci = 0.0

    write(*,*)
    write(*,*) 'Calculating liutex'

    call cpu_time(t5)

    open(18, file='liutexZDZpoints.txt')

    do k = 1, kmax_local
      do j = 1, jmax_local
        do i = 1, imax_local

          ! set velocity gradient tensor
          a(1,1) = dudx(i,j,k)
          a(1,2) = dudy(i,j,k)
          a(1,3) = dudz(i,j,k)
          a(2,1) = dvdx(i,j,k)
          a(2,2) = dvdy(i,j,k)
          a(2,3) = dvdz(i,j,k)
          a(3,1) = dwdx(i,j,k)
          a(3,2) = dwdy(i,j,k)
          a(3,3) = dwdz(i,j,k)

         ! write(*,*) a(3,2)
          !---------------------------------------------------------------------
          ! Cubic Formula
          ! Reference: Numerical Recipes in FORTRAN 77, Second Edition
          ! 5.6 Quadratic and Cubic Equations
          ! Page 179
          !---------------------------------------------------------------------

          ! cubic equation
          ! x**3 + aa * x**2 + bb * x + cc = 0

          ! coefficients of characteristic equation

          aa = -(a(1,1)+a(2,2)+a(3,3))

          tt = matmul(a,a)

          bb = -0.5*(tt(1,1)+tt(2,2)+tt(3,3)-(a(1,1)+a(2,2)+a(3,3))**2)

          cc = -(a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))                          &
               -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))                           &
               +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1)))

          ! discriminant of characteristic equation
          delta = 18*aa*bb*cc-4*aa**3*cc+aa**2*bb**2-4*bb**3-27*cc**2

          qq = (aa**2-3*bb)/9.0
          rr = (2*aa**3-9*aa*bb+27*cc)/54.0

          ! delta = rr**2 - qq**3
          ! alleviate round error
          delta = -delta/108

          vr(1) = 0.0
          vr(2) = 0.0
          vr(3) = 0.0
          if(delta > 0.0) then ! one real root and two complex conjugate roots

            aaaa = -sign(1.0, rr)*(abs(rr)+sqrt(delta))**(1.0/3.0)

            if(aaaa == 0.0) then
              bbbb = 0.0
            else
              bbbb = qq/aaaa
            end if

            eig1c = cmplx(-0.5*(aaaa+bbbb)-aa/3.0, 0.5*sqrt(3.0)*(aaaa-bbbb))
            eig2c = cmplx(real(eig1c), -aimag(eig1c)) !original wrong here
            eig3r = aaaa+bbbb-aa/3.0
            liutexEigr(i,j,k) = abs(eig3r)

            ! real right eigenvector

            delta1 = (a(1,1)-eig3r)*(a(2,2)-eig3r) - a(2,1)*a(1,2)
            delta2 = (a(2,2)-eig3r)*(a(3,3)-eig3r) - a(2,3)*a(3,2)
            delta3 = (a(1,1)-eig3r)*(a(3,3)-eig3r) - a(1,3)*a(3,1)

            if(delta1 == 0.0 .and. delta2 == 0.0 .and. delta3 == 0.0) then
              write(*,*) 'ERROR: delta1 = delta2 = delta3 = 0.0'
              write(*,*) a(1,1)-eig3r,  a(1,2),       a(1,3)
              write(*,*) a(2,1),        a(2,2)-eig3r, a(2,3)
              write(*,*) a(3,1),        a(3,2),       a(3,3)-eig3r
              write(*,*) i, j, k
              stop
            end if

            if(abs(delta1) >= abs(delta2) .and. abs(delta1) >= abs(delta3)) then

              vr(1) = (-(a(2,2)-eig3r)*a(1,3) +         a(1,2)*a(2,3))/delta1
              vr(2) = (         a(2,1)*a(1,3) - (a(1,1)-eig3r)*a(2,3))/delta1
              vr(3) = 1.0

            else if(abs(delta2) >= abs(delta1) .and. abs(delta2) >= abs(delta3)) then

              vr(1) = 1.0
              vr(2) = (-(a(3,3)-eig3r)*a(2,1) +         a(2,3)*a(3,1))/delta2
              vr(3) = (         a(3,2)*a(2,1) - (a(2,2)-eig3r)*a(3,1))/delta2

            else if(abs(delta3) >= abs(delta1) .and. abs(delta3) >= abs(delta2)) then

              vr(1) = (-(a(3,3)-eig3r)*a(1,2) +         a(1,3)*a(3,2))/delta3
              vr(2) = 1.0
              vr(3) = (         a(3,1)*a(1,2) - (a(1,1)-eig3r)*a(3,2))/delta3

            else

              write(*,*) 'ERROR: '
              write(*,*) delta1, delta2, delta3
              stop

            end if

            temp = sqrt(vr(1)**2+vr(2)**2+vr(3)**2)

            vr(1) = vr(1)/temp
            vr(2) = vr(2)/temp
            vr(3) = vr(3)/temp

            call rotation(z0, vr, qqq)

            vg = matmul(transpose(qqq), a)
            vg = matmul(vg, qqq)

            alpha = 0.5*sqrt((vg(2,2)-vg(1,1))**2+(vg(2,1)+vg(1,2))**2)
            beta  = 0.5*(vg(2,1)-vg(1,2))

            if(beta**2 > alpha**2) then

              if(beta > 0.0) then
                rm = 2*(beta-alpha)
                liutex_x(i, j, k) = rm*vr(1)
                liutex_y(i, j, k) = rm*vr(2)
                liutex_z(i, j, k) = rm*vr(3)
              else
                rm = 2*(beta+alpha)
                liutex_x(i, j, k) = rm*vr(1)
                liutex_y(i, j, k) = rm*vr(2)
                liutex_z(i, j, k) = rm*vr(3)
              end if

            else

              liutex_x(i,j,k) = 0.0
              liutex_y(i,j,k) = 0.0
              liutex_z(i,j,k) = 0.0

            end if

            liutex_mag(i,j,k) = sqrt( liutex_x(i,j,k)**2                        &
                                     +liutex_y(i,j,k)**2                        &
                                     +liutex_z(i,j,k)**2 )

            lamb_ci(i,j,k) = abs(aimag(eig1c))

          else ! three real roots

            liutex_x(i,j,k) = 0.0
            liutex_y(i,j,k) = 0.0
            liutex_z(i,j,k) = 0.0
            liutex_mag(i,j,k) = 0.0

          end if

          ljmtmp2 = sqrt( liutex_mag_x(i,j,k)**2 + liutex_mag_y(i,j,k)**2 +liutex_mag_z(i,j,k)**2 )

          isLocExtrOmg(i,j,k) = 0.0

          if (ljmtmp2 < 1e-12 .and. liutex_mag(i,j,k)> 1.0e-5) then
            isLocExtrOmg(i,j,k)=1.0
          end if

          if (isLocExtrOmg(i,j,k) .eq. 1.0) then

            write(18,"(3(f21.12))") x(i+istart-1, j+jstart-1, k+kstart-1),  &
                                    y(i+istart-1, j+jstart-1, k+kstart-1),  &
                                    z(i+istart-1, j+jstart-1, k+kstart-1)

          end if

          ljmtmp = sqrt( (liutex_mag_y(i,j,k)*liutex_z(i,j,k)-liutex_mag_z(i,j,k)*liutex_y(i,j,k))**2 &
                       +(liutex_mag_z(i,j,k)*liutex_x(i,j,k)-liutex_mag_x(i,j,k)*liutex_z(i,j,k))**2 &
                       +(liutex_mag_x(i,j,k)*liutex_y(i,j,k)-liutex_mag_y(i,j,k)*liutex_x(i,j,k))**2)

          if (ljmtmp2>1.0e-6) then
            liutexEigr(i,j,k) = abs(liutex_mag_x(i,j,k)*vr(1)+liutex_mag_y(i,j,k)*vr(2)+liutex_mag_z(i,j,k)*vr(3) )/ljmtmp2
          else
            liutexEigr(i,j,k) = 0.0
          end if

          omgGradXliutex(i,j,k)=1.0

          if(omega(i,j,k)<0.51) then
            omgGradXliutex(i,j,k) = 0.0
          endif
          
        end do
      end do
    end do

    call cpu_time(t6)

    write(*,*) 'calculation time: ', t6 - t5

    !---------------------------------------------------------------------------
    ! Calculating Omega and Modified-Omega-Liutex
    !---------------------------------------------------------------------------

    allocate(omega_l(imax_local, jmax_local, kmax_local))

    allocate(o_alpha(imax_local, jmax_local, kmax_local))
    allocate(o_beta(imax_local, jmax_local, kmax_local))

    allocate(lambda_cr(imax_local, jmax_local, kmax_local))
    allocate(lambda_r(imax_local, jmax_local, kmax_local))

    write(*,*)
    write(*,*) 'Calculating Omega and Modified-Omega-Liutex'

    call cpu_time(t5)

    do k = 1, kmax_local
      do j = 1, jmax_local
        do i = 1, imax_local

          a(1,1) = dudx(i,j,k)
          a(1,2) = dudy(i,j,k)
          a(1,3) = dudz(i,j,k)
          a(2,1) = dvdx(i,j,k)
          a(2,2) = dvdy(i,j,k)
          a(2,3) = dvdz(i,j,k)
          a(3,1) = dwdx(i,j,k)
          a(3,2) = dwdy(i,j,k)
          a(3,3) = dwdz(i,j,k)

          !---------------------------------------------------------------------
          ! Cubic Formula
          ! Reference: Numerical Recipes in FORTRAN 77, Second Edition
          ! 5.6 Quadratic and Cubic Equations
          ! Page 179
          !---------------------------------------------------------------------

          ! cubic equation
          ! x**3 + aa * x**2 + bb * x + cc = 0

          ! coefficients of characteristic equation

          aa = -(a(1,1)+a(2,2)+a(3,3))

          tt = matmul(a,a)

          bb = -0.5*(tt(1,1)+tt(2,2)+tt(3,3)-(a(1,1)+a(2,2)+a(3,3))**2)

          cc = -(a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))                          &
                 -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))                         &
                 +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1)))

          ! discriminant of characteristic equation
          delta = 18*aa*bb*cc-4*aa**3*cc+aa**2*bb**2-4*bb**3-27*cc**2

          qq = (aa**2-3*bb)/9.0
          rr = (2*aa**3-9*aa*bb+27*cc)/54.0

          ! delta = rr**2 - qq**3
          ! alleviate round error
          delta = -delta/108

          if(delta > 0.0) then ! one real root and two complex conjugate roots

            aaaa = -sign(1.0, rr)*(abs(rr)+sqrt(delta))**(1.0/3.0)

            if(aaaa == 0.0) then
              bbbb = 0.0
            else
              bbbb = qq/aaaa
            end if

            eig1c = cmplx(-0.5*(aaaa+bbbb)-aa/3.0, 0.5*sqrt(3.0)*(aaaa-bbbb))
            eig2c = cmplx(real(eig1c), -aimag(eig1c))
            eig3r = aaaa+bbbb-aa/3.0

            ! real right eigenvector

            delta1 = (a(1,1)-eig3r)*(a(2,2)-eig3r) - a(2,1)*a(1,2)
            delta2 = (a(2,2)-eig3r)*(a(3,3)-eig3r) - a(2,3)*a(3,2)
            delta3 = (a(1,1)-eig3r)*(a(3,3)-eig3r) - a(1,3)*a(3,1)

            if(delta1 == 0.0 .and. delta2 == 0.0 .and. delta3 == 0.0) then
              write(*,*) 'ERROR: delta1 = delta2 = delta3 = 0.0'
              write(*,*) a(1,1)-eig3r,  a(1,2),       a(1,3)
              write(*,*) a(2,1),        a(2,2)-eig3r, a(2,3)
              write(*,*) a(3,1),        a(3,2),       a(3,3)-eig3r
              write(*,*) i, j, k
              stop
            end if

            if(abs(delta1) >= abs(delta2) .and.                                &
               abs(delta1) >= abs(delta3)) then

              vr(1) = (-(a(2,2)-eig3r)*a(1,3) +         a(1,2)*a(2,3))/delta1
              vr(2) = (         a(2,1)*a(1,3) - (a(1,1)-eig3r)*a(2,3))/delta1
              vr(3) = 1.0

            else if(abs(delta2) >= abs(delta1) .and.                           &
                    abs(delta2) >= abs(delta3)) then

              vr(1) = 1.0
              vr(2) = (-(a(3,3)-eig3r)*a(2,1) +         a(2,3)*a(3,1))/delta2
              vr(3) = (         a(3,2)*a(2,1) - (a(2,2)-eig3r)*a(3,1))/delta2

            else if(abs(delta3) >= abs(delta1) .and.                           &
                    abs(delta3) >= abs(delta2)) then

              vr(1) = (-(a(3,3)-eig3r)*a(1,2) +         a(1,3)*a(3,2))/delta3
              vr(2) = 1.0
              vr(3) = (         a(3,1)*a(1,2) - (a(1,1)-eig3r)*a(3,2))/delta3

            else

              write(*,*) 'ERROR: '
              write(*,*) delta1, delta2, delta3
              stop

            end if

            temp = sqrt(vr(1)**2+vr(2)**2+vr(3)**2)

            vr(1) = vr(1)/temp
            vr(2) = vr(2)/temp
            vr(3) = vr(3)/temp

            call rotation(z0, vr, qqq)

            vg = matmul(transpose(qqq), a)
            vg = matmul(vg, qqq)

            o_alpha(i,j,k) = 0.5*sqrt((vg(2,2)-vg(1,1))**2+(vg(2,1)+vg(1,2))**2)
            o_beta(i,j,k)  = 0.5*(vg(2,1)-vg(1,2))

            lambda_cr(i,j,k) = real(eig1c)
            lambda_r(i,j,k) = eig3r

          else

            o_alpha(i,j,k) = 0.0
            o_beta(i,j,k)  = 0.0

          end if

          maxbeta_alpha = max(maxbeta_alpha, o_beta(i,j,k)**2-o_alpha(i,j,k)**2)

        end do
      end do
    end do

    do k = 1, kmax_local
      do j = 1, jmax_local
        do i = 1, imax_local

          omega_l(i ,j, k) = o_beta(i,j,k)**2/(o_beta(i,j,k)**2+o_alpha(i,j,k)**2+   &
                             lambda_cr(i,j,k)**2+0.5*lambda_r(i,j,k)**2+       &
                             omega_l_eps*maxbeta_alpha)

        end do
      end do
    end do

    call cpu_time(t6)

    write(*,*) 'Finished Modified Liutex Omega'
    write(*,*) 'calculation time: ', t6 - t5

    ! Writing and Saving Data
    !---------------------------------------------------------------------------
    ! f(1): u velocity
    ! f(2): v velocity
    ! f(3): w velocity
    ! f(4): vorticity x
    ! f(5): vorticity y
    ! f(6): vorticity z
    ! f(7): Omega
    ! f(8): Liutex x
    ! f(9): Liutex y
    ! f(10): Liutex z
    ! f(11): Liutex magnitude
    ! f(12): Q Method
    ! f(13): isLocExtrOmg
    ! f(14): Liutex mag gradient x
    ! f(15): Liutex mag gradient y
    ! f(16): Liutex mag gradient z
    ! f(17): Modified-Omega-Liutex
    !---------------------------------------------------------------------------

    outputfilename1 = "data/"//trim(outputfilename)//'_'//chars//".xyz"
    outputfilename2 = "data/"//trim(outputfilename)//'_'//chars//".fun"

    ! write coordinates
    open(fin1, file=trim(outputfilename1), form='unformatted', action='write')
    write(fin1) imax_local, jmax_local, kmax_local
    write(fin1)                                                                &
               (((x_global(i, j, k), i=istart,iend), j=jstart,jend), k=kstart,kend),  &
               (((y_global(i, j, k), i=istart,iend), j=jstart,jend), k=kstart,kend),  &
               (((z_global(i, j, k), i=istart,iend), j=jstart,jend), k=kstart,kend)

    close(fin1)

    ! Write data

    nvar = 17

    allocate(f(imax_local, jmax_local, kmax_local, nvar))

    f(:,:,:,1) = u
    f(:,:,:,2) = v
    f(:,:,:,3) = w
    f(:,:,:,4) = vorticity_x
    f(:,:,:,5) = vorticity_y
    f(:,:,:,6) = vorticity_z
    f(:,:,:,7) = omega
    f(:,:,:,8) = liutex_x
    f(:,:,:,9) = liutex_y
    f(:,:,:,10) = liutex_z
    f(:,:,:,11) = liutex_mag
    f(:,:,:,12) = q
    f(:,:,:,13) = isLocExtrOmg
    f(:,:,:,14) = liutex_mag_x
    f(:,:,:,15) = liutex_mag_y
    f(:,:,:,16) = liutex_mag_z
    f(:,:,:,17) = omega_l


    open(fin1, file=trim(outputfilename2), form='unformatted', action='write')

    write(fin1) imax_local, jmax_local, kmax_local, nvar

    write(fin1) ((((f(i, j, k, nx), i=1,imax_local), &
                                    j=1,jmax_local), &
                                    k=1,kmax_local), &
                                    nx=1,nvar)

    close(fin1)
    close(18)

    deallocate(f)

    ! Creating Name file for Tecplot
    ! This gives each variable a name in tecplot instead of the default placeholder.

    outputfilename3 = "data/"//trim(outputfilename)//'_'//chars//".nam"

    open(fin1, file=trim(outputfilename3), form='formatted', action='write')

    write(fin1,*) "U"
    write(fin1,*) "V"
    write(fin1,*) "W"
    write(fin1,*) "Vorticity X"
    write(fin1,*) "Vorticity Y"
    write(fin1,*) "Vorticity Z"
    write(fin1,*) "Omega"
    write(fin1,*) "Liutex X"
    write(fin1,*) "Liutex Y"
    write(fin1,*) "Liutex Z"
    write(fin1,*) "Liutex Magnitude"
    write(fin1,*) "Q"
    write(fin1,*) "isLocExtrOmg"
    write(fin1,*) "Liutex Magnitude Gradient X"
    write(fin1,*) "Liutex Magnitude Gradient Y"
    write(fin1,*) "Liutex Magnitude Gradient Z"
    write(fin1,*) "Omega Liutex"

    close(fin1)

    write(*,*) trim(outputfilename)//' finished successfully!'

    deallocate(dudx)
    deallocate(dudy)
    deallocate(dudz)
    deallocate(dvdx)
    deallocate(dvdy)
    deallocate(dvdz)
    deallocate(dwdx)
    deallocate(dwdy)
    deallocate(dwdz)

    deallocate(vorticity_x)
    deallocate(vorticity_y)
    deallocate(vorticity_z)

    deallocate(vorticity_mag)

    deallocate(u)
    deallocate(v)
    deallocate(w)

    deallocate(liutex_x)
    deallocate(liutex_y)
    deallocate(liutex_z)

    deallocate(q)
    deallocate(lamb_ci)
    deallocate(omega)

    deallocate(q_x)
    deallocate(q_y)
    deallocate(q_z)
    deallocate(omgGradXliutex)
    deallocate(isLocExtrOmg)
    deallocate(liutexEigr)
    deallocate(liutex_mag_x)
    deallocate(liutex_mag_y)
    deallocate(liutex_mag_z)
    deallocate(liutex_mag)

    deallocate(o_alpha)
    deallocate(o_beta)

    deallocate(lambda_cr)
    deallocate(lambda_r)

    deallocate(omega_l)

    end do

end program DNSUTA_extract



! calculate the liutex
! a velocity gradient tensor, vor: vorticity,vr: rotational axis, rorMag:rotational strength
subroutine cal_liutex(a, vor, vr,rorMag)
  implicit none

  real, intent(in) :: a(3,3)
  real, intent(in) :: vor(3)
  real, intent(out) :: vr(3), rorMag

  real:: aa,bb,cc,delta,rr,aaaa,bbbb,qq,delta1,delta2,delta3,temp
  complex:: eig1c,eig2c,eig3r
  real :: tt(3,3)
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

          bb = -0.5 * ( tt(1,1) + tt(2,2) + tt(3,3) - ( a(1,1) + a(2,2) + a(3,3) )**2 )

          cc = -(a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))                          &
               -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))                           &
               +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1)))

          ! discriminant of characteristic equation
          delta = 18 * aa * bb * cc - 4 * aa**3 * cc + aa**2 * bb**2 - 4 * bb**3 - 27 * cc**2

          qq = ( aa**2 - 3 * bb ) / 9.0
          rr = ( 2 * aa**3 - 9 * aa * bb + 27 * cc ) / 54.0

          ! delta = rr**2 - qq**3
          ! alleviate round error
          delta = -delta / 108
          vr = 0.0
          rorMag = 0.0

          if(delta > 0.0) then ! one real root and two complex conjugate roots

            aaaa = -sign(1.0, rr) * ( abs(rr) + sqrt(delta) )**(1.0/3.0)

            if(aaaa == 0.0) then
              bbbb = 0.0
            else
              bbbb = qq/aaaa
            end if

            eig1c = cmplx(-0.5 * (aaaa+bbbb) - aa / 3.0, 0.5 * sqrt(3.0) * (aaaa-bbbb))
            eig2c = cmplx(real(eig1c), -aimag(eig1c)) !original wrong here
            eig3r = aaaa + bbbb - aa / 3.0

            ! real right eigenvector

            delta1 = (a(1,1)-eig3r)*(a(2,2)-eig3r) - a(2,1)*a(1,2)
            delta2 = (a(2,2)-eig3r)*(a(3,3)-eig3r) - a(2,3)*a(3,2)
            delta3 = (a(1,1)-eig3r)*(a(3,3)-eig3r) - a(1,3)*a(3,1)

            if(delta1 == 0.0 .and. delta2 == 0.0 .and. delta3 == 0.0) then
              write(*,*) 'ERROR: delta1 = delta2 = delta3 = 0.0'
              write(*,*) a(1,1)-eig3r,  a(1,2),       a(1,3)
              write(*,*) a(2,1),        a(2,2)-eig3r, a(2,3)
              write(*,*) a(3,1),        a(3,2),       a(3,3)-eig3r
            !  write(*,*) i, j, k
              stop
            end if

            if(abs(delta1) >= abs(delta2) .and.                                &
               abs(delta1) >= abs(delta3)) then

              vr(1) = (-(a(2,2)-eig3r)*a(1,3) +         a(1,2)*a(2,3))/delta1
              vr(2) = (         a(2,1)*a(1,3) - (a(1,1)-eig3r)*a(2,3))/delta1
              vr(3) = 1.0

            else if(abs(delta2) >= abs(delta1) .and.                           &
                    abs(delta2) >= abs(delta3)) then

              vr(1) = 1.0
              vr(2) = (-(a(3,3)-eig3r)*a(2,1) +         a(2,3)*a(3,1))/delta2
              vr(3) = (         a(3,2)*a(2,1) - (a(2,2)-eig3r)*a(3,1))/delta2

            else if(abs(delta3) >= abs(delta1) .and.                           &
                    abs(delta3) >= abs(delta2)) then

               vr(1) = (-(a(3,3)-eig3r)*a(1,2) +         a(1,3)*a(3,2))/delta3
               vr(2) = 1.0
               vr(3) = (         a(3,1)*a(1,2) - (a(1,1)-eig3r)*a(3,2))/delta3

            else

              write(*,*) 'ERROR: '
              write(*,*) delta1, delta2, delta3
              stop

            end if

            temp = sqrt(vr(1)**2+vr(2)**2+vr(3)**2)

            vr(1) = vr(1)/temp
            vr(2) = vr(2)/temp
            vr(3) = vr(3)/temp
            temp=dot_product(vor,vr)
            vr(1)=sign(1.0,temp)*vr(1)
            vr(2)=sign(1.0,temp)*vr(2)
            vr(3)=sign(1.0,temp)*vr(3)

            rorMag=(temp - sqrt(temp**2-4*aimag(eig2c)*aimag(eig2c)))
        end if

end subroutine cal_liutex


subroutine rotation(u, v, r)
!-------------------------------------------------------------------------------
! calculate rotation matrix r which rotates unit vector u to unit vector v
!-------------------------------------------------------------------------------

  implicit none

  real, intent(in) :: u(3)
  real, intent(in) :: v(3)
  real, intent(out) :: r(3, 3)

  real :: a(3)
  real :: aa
  real :: t
  real :: alpha
  real :: c, s

  real, parameter :: eps = 1.0e-10

  ! a = u x v
  a(1) = u(2)*v(3)-u(3)*v(2)
  a(2) = u(3)*v(1)-u(1)*v(3)
  a(3) = u(1)*v(2)-u(2)*v(1)

  ! norm
  aa = sqrt(a(1)**2+a(2)**2+a(3)**2)

  if(aa < eps) then

    r = 0.0
    r(1,1) = 1.0
    r(2,2) = 1.0
    r(3,3) = 1.0

  else

    a = a/aa
    t = u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
    if(t > 1.0) t = 1.0
    if(t < -1.0) t = -1.0
    alpha = acos(t)

    c = cos(alpha)
    s = sin(alpha)

    r(1,1) = a(1)**2*(1-c)+c
    r(1,2) = a(1)*a(2)*(1-c)-a(3)*s
    r(1,3) = a(1)*a(3)*(1-c)+a(2)*s

    r(2,1) = a(2)*a(1)*(1-c)+a(3)*s
    r(2,2) = a(2)**2*(1-c)+c
    r(2,3) = a(2)*a(3)*(1-c)-a(1)*s

    r(3,1) = a(3)*a(1)*(1-c)-a(2)*s
    r(3,2) = a(3)*a(2)*(1-c)+a(1)*s
    r(3,3) = a(3)**2*(1-c)+c

  end if

end subroutine rotation


subroutine itoa6(intnum, chars)

  implicit none

  integer :: intnum
  character :: chd1, chd2, chd3, chd4, chd5, chd6
  character(6) :: chars
  integer :: id1, id2, id3, id4, id5, id6

  id1 = (intnum/1     *1      - intnum/10     *10)/1
  id2 = (intnum/10    *10     - intnum/100    *100)/10
  id3 = (intnum/100   *100    - intnum/1000   *1000)/100
  id4 = (intnum/1000  *1000   - intnum/10000  *10000)/1000
  id5 = (intnum/10000 *10000  - intnum/100000 *100000)/10000
  id6 = (intnum/100000*100000 - intnum/1000000*1000000)/100000

  chd1 = char(id1+48)
  chd2 = char(id2+48)
  chd3 = char(id3+48)
  chd4 = char(id4+48)
  chd5 = char(id5+48)
  chd6 = char(id6+48)

  chars = chd6//chd5//chd4//chd3//chd2//chd1

end subroutine



