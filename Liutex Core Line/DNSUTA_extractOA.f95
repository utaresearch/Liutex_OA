program DNSUTA_extract
!-------------------------------------------------------------------------------
! Email: gaoyisheng@me.com, cliu@uta.edu
! Email: oscar.alvarez@mavs.uta.edu
! Department of Mathematics, University of Texas at Arlington,
! Arlington, Texas, USA
! Fall 2021
!-------------------------------------------------------------------------------
  implicit none

  integer :: imax, jmax, kmax
  integer :: imax_local, jmax_local, kmax_local

  integer :: istart, iend
  integer :: jstart, jend
  integer :: kstart, kend

  real(8), dimension(:,:,:), allocatable :: x
  real(8), dimension(:,:,:), allocatable :: y
  real(8), dimension(:,:,:), allocatable :: z

  real(8), dimension(:,:,:), allocatable :: u
  real(8), dimension(:,:,:), allocatable :: v
  real(8), dimension(:,:,:), allocatable :: w

  real(8), dimension(:,:,:,:), allocatable :: f

  real(8), dimension(:,:,:), allocatable :: dudx, dudy, dudz
  real(8), dimension(:,:,:), allocatable :: dvdx, dvdy, dvdz
  real(8), dimension(:,:,:), allocatable :: dwdx, dwdy, dwdz
  real(8) :: rdudx,rdudy,rdudz,rdvdx,rdvdy,rdvdz,rdwdx
  real(8) :: rdwdy,rdwdz,rvort_x,rvort_y,rvort_z

  real(8), dimension(:,:,:), allocatable :: vorticity_x
  real(8), dimension(:,:,:), allocatable :: vorticity_y
  real(8), dimension(:,:,:), allocatable :: vorticity_z

  real(8), dimension(:,:,:), allocatable :: vorticity_mag

  real(8), dimension(:,:,:), allocatable :: omega,localOmega,Qmethod,QGlob
  real(8), dimension(:,:,:), allocatable :: omg_x,omg_y,omg_z, omgGradXRortex
  real(8), dimension(:,:,:), allocatable :: omgLoc_x,omgLoc_y,omgLoc_z,QLoc_x,QLoc_y,QLoc_z
  real(8), dimension(:,:,:), allocatable :: omg_xx,omg_xy,omg_xz
  real(8), dimension(:,:,:), allocatable :: omg_yx,omg_yy,omg_yz
  real(8), dimension(:,:,:), allocatable :: omg_zx,omg_zy,omg_zz

  real(8), dimension(:,:,:), allocatable :: rortex_x
  real(8), dimension(:,:,:), allocatable :: rortex_y
  real(8), dimension(:,:,:), allocatable :: rortex_z
  real(8), dimension(:,:,:), allocatable :: isLocExtrOmg
  real(8), dimension(:,:,:), allocatable :: rortexEigr

  real(8), dimension(:,:,:), allocatable :: rortex_mag, lamb_ci,rortexGlob_mag
  real(8),dimension (:,:,:), allocatable :: rortex_mag_x,rortex_mag_y, rortex_mag_z

  real(8) :: vor(3),vrtmp(3), rorMag,normljm,normljm2,ljmtmp,ljmtmp2

  character(100) :: gridfilename, funcfilename, outputfilename
  character(100) :: outputfilename1
  character(100) :: outputfilename2
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

  real(8) :: aaa
  real(8) :: bbb

  real(8) :: thljm, qljm,rljm,pi
  real(8) :: omega_eps

  real(8) :: vort_2

  real(8) :: u_xi, u_eta, u_zeta, u_xiR, u_etaR,u_zetaR
  real(8) :: v_xi, v_eta, v_zeta, u_xiQ, u_etaQ,u_zetaQ
  real(8) :: w_xi, w_eta, w_zeta

  real(8) :: x_xi, x_eta, x_zeta
  real(8) :: y_xi, y_eta, y_zeta
  real(8) :: z_xi, z_eta, z_zeta

  real(8) :: xi_x, xi_y, xi_z
  real(8) :: eta_x, eta_y, eta_z
  real(8) :: zeta_x, zeta_y, zeta_z

  real(8) :: det

  real(8) :: a(3,3), Aljm(3,3), Matljm(3,3)

  real(8) :: t1, t2, t3, t4, t5, t6

  real(8) :: aa, bb, cc
  real(8) :: delta
  real(8) :: tt(3, 3)

  complex(8) :: eig1c, eig2c
  real(8) :: eig3r

  real(8) :: qq, rr, root11,root22,root33
  real(8) :: aaaa, bbbb


  real(8) :: vr(3)
  real(8) :: temp
  real(8) :: vg(3,3)

  real(8) :: alpha, beta

  ! Modified Omega Liutex variables
  real(8), dimension(:,:,:), allocatable :: omega_l
  real(8), dimension(:,:,:), allocatable :: lambda_cr, lambda_r
  real(8) :: maxbeta_alpha, omega_l_eps

  real(8), dimension(:,:,:), allocatable :: o_alpha, o_beta

  real(8), parameter :: z0(3) = (/0.0, 0.0, 1.0/)
  real(8) :: rm
  real(8) :: qqq(3,3)

  real(8) :: delta1, delta2, delta3

  character(6) :: chars

  continue

  ! the default epsilon value for Omega method
  omega_eps = 1.0e-4
  pi = 4.0*atan(1.0)

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

  allocate(x(imax, jmax, kmax))
  allocate(y(imax, jmax, kmax))
  allocate(z(imax, jmax, kmax))

  write(*,*)
  write(*,*) 'imax = ', imax, ' jmax = ', jmax, ' kmax = ', kmax

  ! read coordinates
  read(fin1, iostat=ios, iomsg=msg)                                            &
                              (((x(i, j, k), i=1,imax), j=1,jmax), k=1,kmax),  &
                              (((y(i, j, k), i=1,imax), j=1,jmax), k=1,kmax),  &
                              (((z(i, j, k), i=1,imax), j=1,jmax), k=1,kmax)

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

    allocate(localOmega(imax_local, jmax_local, kmax_local))
    allocate(omega(imax, jmax, kmax))

    allocate(QGlob(imax, jmax, kmax))
    allocate(Qmethod(imax_local, jmax_local, kmax_local))
    allocate(QLoc_x(imax_local, jmax_local, kmax_local))
    allocate(QLoc_y(imax_local, jmax_local, kmax_local))
    allocate(QLoc_z(imax_local, jmax_local, kmax_local))

    allocate(omg_x(imax, jmax, kmax))

    allocate(omgGradXRortex(imax_local, jmax_local, kmax_local))
    allocate(omgLoc_x(imax_local, jmax_local, kmax_local))
    allocate(omgLoc_y(imax_local, jmax_local, kmax_local))
    allocate(omgLoc_z(imax_local, jmax_local, kmax_local))

    allocate(omg_xx(imax_local, jmax_local, kmax_local))
    allocate(omg_xy(imax_local, jmax_local, kmax_local))
    allocate(omg_xz(imax_local, jmax_local, kmax_local))

    allocate(omg_y(imax, jmax, kmax))
    allocate(omg_yx(imax_local, jmax_local, kmax_local))
    allocate(omg_yy(imax_local, jmax_local, kmax_local))
    allocate(omg_yz(imax_local, jmax_local, kmax_local))

    allocate(omg_z(imax, jmax, kmax))
    allocate(omg_zx(imax_local, jmax_local, kmax_local))
    allocate(omg_zy(imax_local, jmax_local, kmax_local))
    allocate(omg_zz(imax_local, jmax_local, kmax_local))

    allocate(isLocExtrOmg(imax_local, jmax_local, kmax_local))
    allocate(u(imax_local, jmax_local, kmax_local))
    allocate(v(imax_local, jmax_local, kmax_local))
    allocate(w(imax_local, jmax_local, kmax_local))

    allocate(rortex_x(imax_local, jmax_local, kmax_local))
    allocate(rortex_y(imax_local, jmax_local, kmax_local))
    allocate(rortex_z(imax_local, jmax_local, kmax_local))
    allocate(rortexEigr(imax_local, jmax_local, kmax_local))

    allocate(rortex_mag(imax_local, jmax_local, kmax_local))
    allocate(rortex_mag_x(imax_local, jmax_local, kmax_local))
    allocate(rortex_mag_y(imax_local, jmax_local, kmax_local))
    allocate(rortex_mag_z(imax_local, jmax_local, kmax_local))
    allocate(lamb_ci(imax_local, jmax_local, kmax_local))
    allocate(rortexGlob_mag(imax, jmax, kmax))

    isLocExtrOmg = 0.0
    omgGradXRortex = 0.0
    rortexEigr = 0.0

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
    localOmega = 0.0
    Qmethod = 0.0

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

    do k = kstart, kend
      do j = jstart, jend
        do i = istart, iend

          if(i == 1) then

            ! forward difference

            u_xi = f(2, j, k, 2) - f(1, j, k, 2)
            v_xi = f(2, j, k, 3) - f(1, j, k, 3)
            w_xi = f(2, j, k, 4) - f(1, j, k, 4)

            x_xi = x(2, j, k) - x(1, j, k)
            y_xi = y(2, j, k) - y(1, j, k)
            z_xi = z(2, j, k) - z(1, j, k)

          else if(i == imax) then

            ! backward difference

            u_xi = f(imax, j, k, 2)-f(imax-1, j, k, 2)
            v_xi = f(imax, j, k, 3)-f(imax-1, j, k, 3)
            w_xi = f(imax, j, k, 4)-f(imax-1, j, k, 4)

            x_xi = x(imax, j, k)-x(imax-1, j, k)
            y_xi = y(imax, j, k)-y(imax-1, j, k)
            z_xi = z(imax, j, k)-z(imax-1, j, k)

          else

            ! central difference

            u_xi = 0.5*(f(i+1, j, k, 2) - f(i-1, j, k, 2))
            v_xi = 0.5*(f(i+1, j, k, 3) - f(i-1, j, k, 3))
            w_xi = 0.5*(f(i+1, j, k, 4) - f(i-1, j, k, 4))

            x_xi = 0.5*(x(i+1, j, k) - x(i-1, j, k))
            y_xi = 0.5*(y(i+1, j, k) - y(i-1, j, k))
            z_xi = 0.5*(z(i+1, j, k) - z(i-1, j, k))

          end if

          if(j == 1) then

          ! forward difference

          u_eta = f(i, 2, k, 2) - f(i, 1, k, 2)
          v_eta = f(i, 2, k, 3) - f(i, 1, k, 3)
          w_eta = f(i, 2, k, 4) - f(i, 1, k, 4)

          x_eta = x(i, 2, k) - x(i, 1, k)
          y_eta = y(i, 2, k) - y(i, 1, k)
          z_eta = z(i, 2, k) - z(i, 1, k)

          else if(j == jmax) then

          ! backward difference

          u_eta = f(i, jmax, k, 2) - f(i, jmax-1, k, 2)
          v_eta = f(i, jmax, k, 3) - f(i, jmax-1, k, 3)
          w_eta = f(i, jmax, k, 4) - f(i, jmax-1, k, 4)

          x_eta = x(i, jmax, k) - x(i, jmax-1, k)
          y_eta = y(i, jmax, k) - y(i, jmax-1, k)
          z_eta = z(i, jmax, k) - z(i, jmax-1, k)

          else

          ! central difference

          u_eta = 0.5*(f(i, j+1, k, 2) - f(i, j-1, k, 2))
          v_eta = 0.5*(f(i, j+1, k, 3) - f(i, j-1, k, 3))
          w_eta = 0.5*(f(i, j+1, k, 4) - f(i, j-1, k, 4))

          x_eta = 0.5*(x(i, j+1, k) - x(i, j-1, k))
          y_eta = 0.5*(y(i, j+1, k) - y(i, j-1, k))
          z_eta = 0.5*(z(i, j+1, k) - z(i, j-1, k))

          end if

          if(k == 1) then

          ! forward difference

          u_zeta = f(i, j, 2, 2)-f(i, j, 1, 2)
          v_zeta = f(i, j, 2, 3)-f(i, j, 1, 3)
          w_zeta = f(i, j, 2, 4)-f(i, j, 1, 4)

          x_zeta = x(i, j, 2) - x(i, j, 1)
          y_zeta = y(i, j, 2) - y(i, j, 1)
          z_zeta = z(i, j, 2) - z(i, j, 1)

          else if(k == kmax) then

          ! backward difference

          u_zeta = f(i, j, kmax, 2)-f(i, j, kmax-1, 2)
          v_zeta = f(i, j, kmax, 3)-f(i, j, kmax-1, 3)
          w_zeta = f(i, j, kmax, 4)-f(i, j, kmax-1, 4)

          x_zeta = x(i, j, kmax) - x(i, j, kmax-1)
          y_zeta = y(i, j, kmax) - y(i, j, kmax-1)
          z_zeta = z(i, j, kmax) - z(i, j, kmax-1)

          else

          ! central difference

          u_zeta = 0.5*(f(i, j, k+1, 2) - f(i, j, k-1, 2))
          v_zeta = 0.5*(f(i, j, k+1, 3) - f(i, j, k-1, 3))
          w_zeta = 0.5*(f(i, j, k+1, 4) - f(i, j, k-1, 4))

          x_zeta = 0.5*(x(i, j, k+1) - x(i, j, k-1))
          y_zeta = 0.5*(y(i, j, k+1) - y(i, j, k-1))
          z_zeta = 0.5*(z(i, j, k+1) - z(i, j, k-1))

          end if

          ! determinant of Jacobian
          det =   x_xi*(y_eta*z_zeta-y_zeta*z_eta)                             &
                - x_eta*(y_xi*z_zeta-y_zeta*z_xi)                              &
                + x_zeta*(y_xi*z_eta-y_eta*z_xi)

          det = 1.0/det

          xi_x = det*(y_eta*z_zeta - y_zeta*z_eta)
          xi_y = det*(x_zeta*z_eta - x_eta*z_zeta)
          xi_z = det*(x_eta*y_zeta - x_zeta*y_eta)

          eta_x = det*(y_zeta*z_xi - y_xi*z_zeta)
          eta_y = det*(x_xi*z_zeta - x_zeta*z_xi)
          eta_z = det*(x_zeta*y_xi - x_xi*y_zeta)

          zeta_x = det*(y_xi*z_eta - y_eta*z_xi)
          zeta_y = det*(x_eta*z_xi - x_xi*z_eta)
          zeta_z = det*(x_xi*y_eta - x_eta*y_xi)


          ii = i-istart+1
          jj = j-jstart+1
          kk = k-kstart+1

          ! velocity gradient
          dudx(ii, jj, kk) = u_xi*xi_x + u_eta*eta_x + u_zeta*zeta_x
          dudy(ii, jj, kk) = u_xi*xi_y + u_eta*eta_y + u_zeta*zeta_y
          dudz(ii, jj, kk) = u_xi*xi_z + u_eta*eta_z + u_zeta*zeta_z

          dvdx(ii, jj, kk) = v_xi*xi_x + v_eta*eta_x + v_zeta*zeta_x
          dvdy(ii, jj, kk) = v_xi*xi_y + v_eta*eta_y + v_zeta*zeta_y
          dvdz(ii, jj, kk) = v_xi*xi_z + v_eta*eta_z + v_zeta*zeta_z

          dwdx(ii, jj, kk) = w_xi*xi_x + w_eta*eta_x + w_zeta*zeta_x
          dwdy(ii, jj, kk) = w_xi*xi_y + w_eta*eta_y + w_zeta*zeta_y
          dwdz(ii, jj, kk) = w_xi*xi_z + w_eta*eta_z + w_zeta*zeta_z

          aaa = dudx(ii, jj, kk)**2+dvdy(ii, jj, kk)**2+dwdz(ii, jj, kk)**2 +  &
                0.5*((dudy(ii, jj, kk)+dvdx(ii, jj, kk))**2 +                  &
                   (dudz(ii, jj, kk)+dwdx(ii, jj, kk))**2 +                    &
                   (dvdz(ii, jj, kk)+dwdy(ii, jj, kk))**2)

          ! vorticity
          vorticity_x(ii,jj,kk) = dwdy(ii, jj, kk)-dvdz(ii, jj, kk)
          vorticity_y(ii,jj,kk) = dudz(ii, jj, kk)-dwdx(ii, jj, kk)
          vorticity_z(ii,jj,kk) = dvdx(ii, jj, kk)-dudy(ii, jj, kk)

          vort_2 = vorticity_x(ii,jj,kk)**2                                    &
                  +vorticity_y(ii,jj,kk)**2                                    &
                  +vorticity_z(ii,jj,kk)**2

          vorticity_mag(ii,jj,kk) = sqrt(vort_2)

          bbb = 0.5*vort_2

          localOmega(ii, jj, kk) = bbb/(aaa+bbb+omega_eps)
          Qmethod(ii, jj, kk) = (bbb-aaa)/2.0


          u(ii, jj, kk) = f(i, j, k, 2)
          v(ii, jj, kk) = f(i, j, k, 3)
          w(ii, jj, kk) = f(i ,j, k, 4)

        end do
      end do
    end do

    !--------


    do k = 1, kmax
      do j = 1, jmax
        do i = 1, imax

          if(i == 1) then

            ! forward difference

            u_xi = f(2, j, k, 2) - f(1, j, k, 2)
            v_xi = f(2, j, k, 3) - f(1, j, k, 3)
            w_xi = f(2, j, k, 4) - f(1, j, k, 4)

            x_xi = x(2, j, k) - x(1, j, k)
            y_xi = y(2, j, k) - y(1, j, k)
            z_xi = z(2, j, k) - z(1, j, k)

          else if(i == imax) then

            ! backward difference

            u_xi = f(imax, j, k, 2)-f(imax-1, j, k, 2)
            v_xi = f(imax, j, k, 3)-f(imax-1, j, k, 3)
            w_xi = f(imax, j, k, 4)-f(imax-1, j, k, 4)

            x_xi = x(imax, j, k)-x(imax-1, j, k)
            y_xi = y(imax, j, k)-y(imax-1, j, k)
            z_xi = z(imax, j, k)-z(imax-1, j, k)

          else

            ! central difference

            u_xi = 0.5*(f(i+1, j, k, 2) - f(i-1, j, k, 2))
            v_xi = 0.5*(f(i+1, j, k, 3) - f(i-1, j, k, 3))
            w_xi = 0.5*(f(i+1, j, k, 4) - f(i-1, j, k, 4))

            x_xi = 0.5*(x(i+1, j, k) - x(i-1, j, k))
            y_xi = 0.5*(y(i+1, j, k) - y(i-1, j, k))
            z_xi = 0.5*(z(i+1, j, k) - z(i-1, j, k))

          end if

          if(j == 1) then

            ! forward difference

            u_eta = f(i, 2, k, 2) - f(i, 1, k, 2)
            v_eta = f(i, 2, k, 3) - f(i, 1, k, 3)
            w_eta = f(i, 2, k, 4) - f(i, 1, k, 4)

            x_eta = x(i, 2, k) - x(i, 1, k)
            y_eta = y(i, 2, k) - y(i, 1, k)
            z_eta = z(i, 2, k) - z(i, 1, k)

          else if(j == jmax) then

            ! backward difference

            u_eta = f(i, jmax, k, 2) - f(i, jmax-1, k, 2)
            v_eta = f(i, jmax, k, 3) - f(i, jmax-1, k, 3)
            w_eta = f(i, jmax, k, 4) - f(i, jmax-1, k, 4)

            x_eta = x(i, jmax, k) - x(i, jmax-1, k)
            y_eta = y(i, jmax, k) - y(i, jmax-1, k)
            z_eta = z(i, jmax, k) - z(i, jmax-1, k)

          else

            ! central difference

            u_eta = 0.5*(f(i, j+1, k, 2) - f(i, j-1, k, 2))
            v_eta = 0.5*(f(i, j+1, k, 3) - f(i, j-1, k, 3))
            w_eta = 0.5*(f(i, j+1, k, 4) - f(i, j-1, k, 4))

            x_eta = 0.5*(x(i, j+1, k) - x(i, j-1, k))
            y_eta = 0.5*(y(i, j+1, k) - y(i, j-1, k))
            z_eta = 0.5*(z(i, j+1, k) - z(i, j-1, k))

          end if

          if(k == 1) then

            ! forward difference

            u_zeta = f(i, j, 2, 2)-f(i, j, 1, 2)
            v_zeta = f(i, j, 2, 3)-f(i, j, 1, 3)
            w_zeta = f(i, j, 2, 4)-f(i, j, 1, 4)

            x_zeta = x(i, j, 2) - x(i, j, 1)
            y_zeta = y(i, j, 2) - y(i, j, 1)
            z_zeta = z(i, j, 2) - z(i, j, 1)

          else if(k == kmax) then

            ! backward difference

            u_zeta = f(i, j, kmax, 2)-f(i, j, kmax-1, 2)
            v_zeta = f(i, j, kmax, 3)-f(i, j, kmax-1, 3)
            w_zeta = f(i, j, kmax, 4)-f(i, j, kmax-1, 4)

            x_zeta = x(i, j, kmax) - x(i, j, kmax-1)
            y_zeta = y(i, j, kmax) - y(i, j, kmax-1)
            z_zeta = z(i, j, kmax) - z(i, j, kmax-1)

          else

            ! central difference

            u_zeta = 0.5*(f(i, j, k+1, 2) - f(i, j, k-1, 2))
            v_zeta = 0.5*(f(i, j, k+1, 3) - f(i, j, k-1, 3))
            w_zeta = 0.5*(f(i, j, k+1, 4) - f(i, j, k-1, 4))

            x_zeta = 0.5*(x(i, j, k+1) - x(i, j, k-1))
            y_zeta = 0.5*(y(i, j, k+1) - y(i, j, k-1))
            z_zeta = 0.5*(z(i, j, k+1) - z(i, j, k-1))

          end if

          ! determinant of Jacobian
          det =   x_xi*(y_eta*z_zeta-y_zeta*z_eta)                             &
                - x_eta*(y_xi*z_zeta-y_zeta*z_xi)                              &
                + x_zeta*(y_xi*z_eta-y_eta*z_xi)

          det = 1.0/det

          xi_x = det*(y_eta*z_zeta - y_zeta*z_eta)
          xi_y = det*(x_zeta*z_eta - x_eta*z_zeta)
          xi_z = det*(x_eta*y_zeta - x_zeta*y_eta)

          eta_x = det*(y_zeta*z_xi - y_xi*z_zeta)
          eta_y = det*(x_xi*z_zeta - x_zeta*z_xi)
          eta_z = det*(x_zeta*y_xi - x_xi*y_zeta)

          zeta_x = det*(y_xi*z_eta - y_eta*z_xi)
          zeta_y = det*(x_eta*z_xi - x_xi*z_eta)
          zeta_z = det*(x_xi*y_eta - x_eta*y_xi)

          ! velocity gradient
          rdudx = u_xi*xi_x + u_eta*eta_x + u_zeta*zeta_x
          rdudy = u_xi*xi_y + u_eta*eta_y + u_zeta*zeta_y
          rdudz = u_xi*xi_z + u_eta*eta_z + u_zeta*zeta_z

          rdvdx = v_xi*xi_x + v_eta*eta_x + v_zeta*zeta_x
          rdvdy = v_xi*xi_y + v_eta*eta_y + v_zeta*zeta_y
          rdvdz = v_xi*xi_z + v_eta*eta_z + v_zeta*zeta_z

          rdwdx = w_xi*xi_x + w_eta*eta_x + w_zeta*zeta_x
          rdwdy = w_xi*xi_y + w_eta*eta_y + w_zeta*zeta_y
          rdwdz = w_xi*xi_z + w_eta*eta_z + w_zeta*zeta_z

          aaa = rdudx**2+rdvdy**2+rdwdz**2 +  &
                0.5*((rdudy+rdvdx)**2 +                  &
                   (rdudz+rdwdx)**2 +                    &
                   (rdvdz+rdwdy)**2)

          ! vorticity
          rvort_x = rdwdy-rdvdz
          rvort_y = rdudz-rdwdx
          rvort_z = rdvdx-rdudy

          vor(1)=rvort_x
          vor(2)=rvort_y
          vor(3)=rvort_z

          vort_2 = rvort_x**2                                    &
                  +rvort_y**2                                    &
                  +rvort_z**2

          bbb = 0.5*vort_2

          omega(i, j, k) = bbb/(aaa+bbb+omega_eps)

          a(1,1)=rdudx
          a(1,2)=rdudy
          a(1,3)=rdudz
          a(2,1)=rdvdx
          a(2,2)=rdvdy
          a(2,3)=rdvdz
          a(3,1)=rdwdx
          a(3,2)=rdwdy
          a(3,3)=rdwdz

          call cal_rortex(a, vor, vrtmp,rorMag)

          rortexGlob_mag(i,j,k)=rorMag
          QGlob(i,j,k)=0.5*(bbb-aaa)

        end do
      end do
    end do

    !+++++++

    call cpu_time(t2)

    write(*,*) 'calculation time: ', t2-t1

  ! calculate hessain matrix of the omega
  ! fist step to calculate the first order dirivative of the omega

    do k = 1,kmax
      do j = 1, jmax
        do i = 1, imax

          if(i == 1) then

            ! forward difference

            u_xi = omega(2, j, k) - omega(1, j, k)
            u_xiR =rortexGlob_mag(2, j, k) - rortexGlob_mag(1, j, k)
            u_xiQ =QGlob(2, j, k) - QGlob(1, j, k)

            x_xi = x(2, j, k) - x(1, j, k)
            y_xi = y(2, j, k) - y(1, j, k)
            z_xi = z(2, j, k) - z(1, j, k)

          else if(i == imax) then

            ! backward difference

            u_xi = omega(imax, j, k)-omega(imax-1, j, k)
            u_xiR = rortexGlob_mag(imax, j, k)-rortexGlob_mag(imax-1, j, k)
            u_xiQ = QGlob(imax, j, k)-QGlob(imax-1, j, k)

            x_xi = x(imax, j, k)-x(imax-1, j, k)
            y_xi = y(imax, j, k)-y(imax-1, j, k)
            z_xi = z(imax, j, k)-z(imax-1, j, k)

          else

            ! central difference

            u_xi = 0.5*(omega(i+1, j, k) - omega(i-1, j, k))
            u_xiR = 0.5*(rortexGlob_mag(i+1, j, k) -rortexGlob_mag(i-1, j, k))
            u_xiQ = 0.5*(QGlob(i+1, j, k) -QGlob(i-1, j, k))

            x_xi = 0.5*(x(i+1, j, k) - x(i-1, j, k))
            y_xi = 0.5*(y(i+1, j, k) - y(i-1, j, k))
            z_xi = 0.5*(z(i+1, j, k) - z(i-1, j, k))

          end if

          if(j == 1) then

            ! forward difference

            u_eta = omega(i, 2, k) - omega(i, 1, k)
            u_etaR = rortexGlob_mag(i, 2, k) - rortexGlob_mag(i, 1, k)
            u_etaQ = QGlob(i, 2, k) - QGlob(i, 1, k)

            x_eta = x(i, 2, k) - x(i, 1, k)
            y_eta = y(i, 2, k) - y(i, 1, k)
            z_eta = z(i, 2, k) - z(i, 1, k)

          else if(j == jmax) then

            ! backward difference

            u_eta = omega(i, jmax, k) - omega(i, jmax-1, k)
            u_etaR = rortexGlob_mag(i, jmax, k) - rortexGlob_mag(i, jmax-1, k)
            u_etaQ = QGlob(i, jmax, k) - QGlob(i, jmax-1, k)

            x_eta = x(i, jmax, k) - x(i, jmax-1, k)
            y_eta = y(i, jmax, k) - y(i, jmax-1, k)
            z_eta = z(i, jmax, k) - z(i, jmax-1, k)

          else

            ! central difference

            u_eta = 0.5*(omega(i, j+1, k) - omega(i, j-1, k))
            u_etaR = 0.5*(rortexGlob_mag(i, j+1, k) - rortexGlob_mag(i, j-1, k))
            u_etaQ = 0.5*(QGlob(i, j+1, k) - QGlob(i, j-1, k))

            x_eta = 0.5*(x(i, j+1, k) - x(i, j-1, k))
            y_eta = 0.5*(y(i, j+1, k) - y(i, j-1, k))
            z_eta = 0.5*(z(i, j+1, k) - z(i, j-1, k))

          end if

          if(k == 1) then

            ! forward difference

            u_zeta = omega(i, j, 2)-omega(i, j, 1)
            u_zetaR = rortexGlob_mag(i, j, 2)-rortexGlob_mag(i, j, 1)
            u_zetaQ = QGlob(i, j, 2)-QGlob(i, j, 1)

            x_zeta = x(i, j, 2) - x(i, j, 1)
            y_zeta = y(i, j, 2) - y(i, j, 1)
            z_zeta = z(i, j, 2) - z(i, j, 1)

          else if(k == kmax) then

            ! backward difference

            u_zeta = omega(i, j, kmax)-omega(i, j, kmax-1)
            u_zetaR = rortexGlob_mag(i, j, kmax)-rortexGlob_mag(i, j, kmax-1)
            u_zetaQ = QGlob(i, j, kmax)-QGlob(i, j, kmax-1)

            x_zeta = x(i, j, kmax) - x(i, j, kmax-1)
            y_zeta = y(i, j, kmax) - y(i, j, kmax-1)
            z_zeta = z(i, j, kmax) - z(i, j, kmax-1)

          else

            ! central difference

            u_zeta = 0.5*(omega(i, j, k+1) - omega(i, j, k-1))
            u_zetaR = 0.5*(rortexGlob_mag(i, j, k+1) - rortexGlob_mag(i, j, k-1))
            u_zetaQ = 0.5*(QGlob(i, j, k+1) - QGlob(i, j, k-1))

            x_zeta = 0.5*(x(i, j, k+1) - x(i, j, k-1))
            y_zeta = 0.5*(y(i, j, k+1) - y(i, j, k-1))
            z_zeta = 0.5*(z(i, j, k+1) - z(i, j, k-1))

          end if

          ! determinant of Jacobian
          det =   x_xi*(y_eta*z_zeta-y_zeta*z_eta)                             &
                - x_eta*(y_xi*z_zeta-y_zeta*z_xi)                              &
                + x_zeta*(y_xi*z_eta-y_eta*z_xi)

          det = 1.0/det

          xi_x = det*(y_eta*z_zeta - y_zeta*z_eta)
          xi_y = det*(x_zeta*z_eta - x_eta*z_zeta)
          xi_z = det*(x_eta*y_zeta - x_zeta*y_eta)

          eta_x = det*(y_zeta*z_xi - y_xi*z_zeta)
          eta_y = det*(x_xi*z_zeta - x_zeta*z_xi)
          eta_z = det*(x_zeta*y_xi - x_xi*y_zeta)

          zeta_x = det*(y_xi*z_eta - y_eta*z_xi)
          zeta_y = det*(x_eta*z_xi - x_xi*z_eta)
          zeta_z = det*(x_xi*y_eta - x_eta*y_xi)

          ii = i-istart+1
          jj = j-jstart+1
          kk = k-kstart+1

          ! velocity gradient
          omg_x(i, j, k) = u_xi*xi_x + u_eta*eta_x + u_zeta*zeta_x
          omg_y(i, j, k) = u_xi*xi_y + u_eta*eta_y + u_zeta*zeta_y
          omg_z(i, j, k) = u_xi*xi_z + u_eta*eta_z + u_zeta*zeta_z

          if(i>=istart .and. i <=iend .and. j>=jstart .and. j<= jend .and. k>=kstart .and. k<= kend ) then
              omgLoc_x(ii, jj, kk) = u_xi*xi_x + u_eta*eta_x + u_zeta*zeta_x
              omgLoc_y(ii, jj, kk) = u_xi*xi_y + u_eta*eta_y + u_zeta*zeta_y
              omgLoc_z(ii, jj, kk) = u_xi*xi_z + u_eta*eta_z + u_zeta*zeta_z

              rortex_mag_x(ii, jj, kk) = u_xiR*xi_x + u_etaR*eta_x + u_zetaR*zeta_x
              rortex_mag_y(ii, jj, kk) = u_xiR*xi_y + u_etaR*eta_y + u_zetaR*zeta_y
              rortex_mag_z(ii, jj, kk) = u_xiR*xi_z + u_etaR*eta_z + u_zetaR*zeta_z

              QLoc_x(ii, jj, kk) = u_xiQ*xi_x + u_etaQ*eta_x + u_zetaQ*zeta_x
              QLoc_y(ii, jj, kk) = u_xiQ*xi_y + u_etaQ*eta_y + u_zetaQ*zeta_y
              QLoc_z(ii, jj, kk) = u_xiQ*xi_z + u_etaQ*eta_z + u_zetaQ*zeta_z
          end if

        end do
      end do
    end do

    ! end of the first step of hessain matrix computation

    do k = kstart, kend
      do j = jstart, jend
        do i = istart, iend

          if(i == 1) then

            ! forward difference

            u_xi = omg_x(2, j, k) - omg_x(1, j, k)
            v_xi = omg_y(2, j, k) - omg_y(1, j, k)
            w_xi = omg_z(2, j, k) - omg_z(1, j, k)

            x_xi = x(2, j, k) - x(1, j, k)
            y_xi = y(2, j, k) - y(1, j, k)
            z_xi = z(2, j, k) - z(1, j, k)

          else if(i == imax) then

            ! backward difference

            u_xi = omg_x(imax, j, k)-omg_x(imax-1, j, k)
            v_xi = omg_y(imax, j, k)-omg_y(imax-1, j, k)
            w_xi = omg_z(imax, j, k)-omg_z(imax-1, j, k)

            x_xi = x(imax, j, k)-x(imax-1, j, k)
            y_xi = y(imax, j, k)-y(imax-1, j, k)
            z_xi = z(imax, j, k)-z(imax-1, j, k)

          else

            ! central difference

            u_xi = 0.5d0*(omg_x(i+1, j, k) - omg_x(i-1, j, k))
            v_xi = 0.5d0*(omg_y(i+1, j, k) - omg_y(i-1, j, k))
            w_xi = 0.5d0*(omg_z(i+1, j, k) - omg_z(i-1, j, k))

            x_xi = 0.5d0*(x(i+1, j, k) - x(i-1, j, k))
            y_xi = 0.5d0*(y(i+1, j, k) - y(i-1, j, k))
            z_xi = 0.5d0*(z(i+1, j, k) - z(i-1, j, k))

          end if

          if(j == 1) then

            ! forward difference

            u_eta = omg_x(i, 2, k) - omg_x(i, 1, k)
            v_eta = omg_y(i, 2, k) - omg_y(i, 1, k)
            w_eta = omg_z(i, 2, k) - omg_z(i, 1, k)

            x_eta = x(i, 2, k) - x(i, 1, k)
            y_eta = y(i, 2, k) - y(i, 1, k)
            z_eta = z(i, 2, k) - z(i, 1, k)

          else if(j == jmax) then

            ! backward difference

            u_eta = omg_x(i, jmax, k) - omg_x(i, jmax-1, k)
            v_eta = omg_y(i, jmax, k) - omg_y(i, jmax-1, k)
            w_eta = omg_z(i, jmax, k) - omg_z(i, jmax-1, k)

            x_eta = x(i, jmax, k) - x(i, jmax-1, k)
            y_eta = y(i, jmax, k) - y(i, jmax-1, k)
            z_eta = z(i, jmax, k) - z(i, jmax-1, k)

          else

            ! central difference

            u_eta = 0.5d0*(omg_x(i, j+1, k) - omg_x(i, j-1, k))
            v_eta = 0.5d0*(omg_y(i, j+1, k) - omg_y(i, j-1, k))
            w_eta = 0.5d0*(omg_z(i, j+1, k) - omg_z(i, j-1, k))

            x_eta = 0.5d0*(x(i, j+1, k) - x(i, j-1, k))
            y_eta = 0.5d0*(y(i, j+1, k) - y(i, j-1, k))
            z_eta = 0.5d0*(z(i, j+1, k) - z(i, j-1, k))

          end if

          if(k == 1) then

            ! forward difference

            u_zeta = omg_x(i, j, 2)-omg_x(i, j, 1)
            v_zeta = omg_y(i, j, 2)-omg_y(i, j, 1)
            w_zeta = omg_z(i, j, 2)-omg_z(i, j, 1)

            x_zeta = x(i, j, 2) - x(i, j, 1)
            y_zeta = y(i, j, 2) - y(i, j, 1)
            z_zeta = z(i, j, 2) - z(i, j, 1)

          else if(k == kmax) then

            ! backward difference

            u_zeta = omg_x(i, j, kmax)-omg_x(i, j, kmax-1)
            v_zeta = omg_y(i, j, kmax)-omg_y(i, j, kmax-1)
            w_zeta = omg_z(i, j, kmax)-omg_z(i, j, kmax-1)

            x_zeta = x(i, j, kmax) - x(i, j, kmax-1)
            y_zeta = y(i, j, kmax) - y(i, j, kmax-1)
            z_zeta = z(i, j, kmax) - z(i, j, kmax-1)

          else

            ! central difference

            u_zeta = 0.5d0*(omg_x(i, j, k+1) - omg_x(i, j, k-1))
            v_zeta = 0.5d0*(omg_y(i, j, k+1) - omg_y(i, j, k-1))
            w_zeta = 0.5d0*(omg_z(i, j, k+1) - omg_z(i, j, k-1))

            x_zeta = 0.5d0*(x(i, j, k+1) - x(i, j, k-1))
            y_zeta = 0.5d0*(y(i, j, k+1) - y(i, j, k-1))
            z_zeta = 0.5d0*(z(i, j, k+1) - z(i, j, k-1))

          end if

          ! determinant of Jacobian
          det =   x_xi*(y_eta*z_zeta-y_zeta*z_eta)                             &
                - x_eta*(y_xi*z_zeta-y_zeta*z_xi)                              &
                + x_zeta*(y_xi*z_eta-y_eta*z_xi)

          det = 1.0/det

          xi_x = det*(y_eta*z_zeta - y_zeta*z_eta)
          xi_y = det*(x_zeta*z_eta - x_eta*z_zeta)
          xi_z = det*(x_eta*y_zeta - x_zeta*y_eta)

          eta_x = det*(y_zeta*z_xi - y_xi*z_zeta)
          eta_y = det*(x_xi*z_zeta - x_zeta*z_xi)
          eta_z = det*(x_zeta*y_xi - x_xi*y_zeta)

          zeta_x = det*(y_xi*z_eta - y_eta*z_xi)
          zeta_y = det*(x_eta*z_xi - x_xi*z_eta)
          zeta_z = det*(x_xi*y_eta - x_eta*y_xi)

          ii = i-istart+1
          jj = j-jstart+1
          kk = k-kstart+1

          ! velocity gradient
          omg_xx(ii, jj, kk) = u_xi*xi_x + u_eta*eta_x + u_zeta*zeta_x
          omg_xy(ii, jj, kk) = u_xi*xi_y + u_eta*eta_y + u_zeta*zeta_y
          omg_xz(ii, jj, kk) = u_xi*xi_z + u_eta*eta_z + u_zeta*zeta_z

          omg_yx(ii, jj, kk) = v_xi*xi_x + v_eta*eta_x + v_zeta*zeta_x
          omg_yy(ii, jj, kk) = v_xi*xi_y + v_eta*eta_y + v_zeta*zeta_y
          omg_yz(ii, jj, kk) = v_xi*xi_z + v_eta*eta_z + v_zeta*zeta_z

          omg_zx(ii, jj, kk) = w_xi*xi_x + w_eta*eta_x + w_zeta*zeta_x
          omg_zy(ii, jj, kk) = w_xi*xi_y + w_eta*eta_y + w_zeta*zeta_y
          omg_zz(ii, jj, kk) = w_xi*xi_z + w_eta*eta_z + w_zeta*zeta_z

          !calculate the eigenvalue

          ! set velocity gradient tensor
          a(1,1) = omg_xx(ii, jj, kk)
          a(1,2) = omg_xy(ii, jj, kk)
          a(1,3) = omg_xz(ii, jj, kk)
          a(2,1) = omg_yx(ii, jj, kk)
          a(2,2) = omg_yy(ii, jj, kk)
          a(2,3) = omg_yz(ii, jj, kk)
          a(3,1) = omg_zx(ii, jj, kk)
          a(3,2) = omg_zy(ii, jj, kk)
          a(3,3) = omg_zz(ii, jj, kk)

          Aljm = (a+transpose(a)) / 2.d0

          Matljm = Aljm
          aa = -(Matljm(1,1)+Matljm(2,2)+Matljm(3,3))

          tt = matmul(Matljm,Matljm)

          bb = -0.5d0*(tt(1,1)+tt(2,2)+tt(3,3)-(Matljm(1,1)+Matljm(2,2)+Matljm(3,3))**2)

          cc = -(Matljm(1,1)*(Matljm(2,2)*Matljm(3,3)-Matljm(2,3)*Matljm(3,2))                            &
                 -Matljm(1,2)*(Matljm(2,1)*Matljm(3,3)-Matljm(2,3)*Matljm(3,1))                           &
                 +Matljm(1,3)*(Matljm(2,1)*Matljm(3,2)-Matljm(2,2)*Matljm(3,1)))

          qljm = (aa**2 - 3.d0*bb)/9.d0
          rljm = (2.d0*aa**3 - 9.d0*aa*bb + 27.d0*cc)/54.d0
          thljm = acos(rljm/sqrt(qljm**3))

          root11 = -2.d0*sqrt(qljm)*cos(thljm/3.0d0) - aa/3.d0
          root22 = -2.d0*sqrt(qljm)*cos(thljm/3.0d0 +2.0/3.d0*pi) - aa/3.d0
          root33 = -2.d0*sqrt(qljm)*cos(thljm/3.d0 - 2.0d0/3.d0*pi) - aa/3.d0

          if(root11>root22) then
              qljm=root11
              root11=root22
              root22=qljm
          end if

          if(root22>root33) then
              qljm=root22
              root22=root33
              root33=qljm
          end if

          if(root11>root22) then
              qljm=root11
              root11=root22
              root22=qljm
          end if

        end do
      end do
    end do

    deallocate(omg_x)
    deallocate(omg_y)
    deallocate(omg_z)

    deallocate(omg_xx)
    deallocate(omg_xy)
    deallocate(omg_xz)

    deallocate(omg_yx)
    deallocate(omg_yy)
    deallocate(omg_yz)

    deallocate(omg_zx)
    deallocate(omg_zy)
    deallocate(omg_zz)

    ! end of the hessain matrix computation


    deallocate(f)

    !---------------------------------------------------------------------------
    ! calculate Rortex
    !---------------------------------------------------------------------------

    lamb_ci = 0.0

    write(*,*)
    write(*,*) 'Calculating Rortex'

    call cpu_time(t5)

    do k = 1, kmax_local
      do j = 1, jmax_local
        do i = 1, imax_local

          ! set velocity gradient tensor
          a(1,1) = dudx(i, j, k)
          a(1,2) = dudy(i, j, k)
          a(1,3) = dudz(i, j, k)
          a(2,1) = dvdx(i, j, k)
          a(2,2) = dvdy(i, j, k)
          a(2,3) = dvdz(i, j, k)
          a(3,1) = dwdx(i, j, k)
          a(3,2) = dwdy(i, j, k)
          a(3,3) = dwdz(i, j, k)

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

          bb = -0.5d0*(tt(1,1)+tt(2,2)+tt(3,3)-(a(1,1)+a(2,2)+a(3,3))**2)

          cc = -(a(1,1) * (a(2,2)*a(3,3)-a(2,3)*a(3,2))                          &
               -a(1,2) * (a(2,1)*a(3,3)-a(2,3)*a(3,1))                           &
               +a(1,3) * (a(2,1)*a(3,2)-a(2,2)*a(3,1)))

          ! discriminant of characteristic equation
          delta = 18*aa*bb*cc-4*aa**3*cc+aa**2*bb**2-4*bb**3-27*cc**2

          qq = (aa**2 - 3.d0*bb) / 9.0d0
          rr = (2.d0 * aa**3 - 9.d0 * aa * bb + 27.d0 * cc) / 54.0d0

          ! delta = rr**2 - qq**3
          ! alleviate round error
          delta = -delta / 108.0d0

          vr(1) = 0.0d0
          vr(2) = 0.0d0
          vr(3) = 0.0d0
          
          if(delta > 0.0d0) then ! one real root and two complex conjugate roots

            aaaa = -sign(1.0, rr)*(abs(rr)+sqrt(delta))**(1.0d0/3.0d0)

            if(aaaa == 0.0d0) then
              bbbb = 0.0d0
            else
              bbbb = qq/aaaa
            end if

            eig1c = cmplx(-0.5d0*(aaaa+bbbb)-aa/3.0d0, 0.5d0*sqrt(3.0d0)*(aaaa-bbbb), kind=8)
            eig2c = cmplx(real(eig1c), -aimag(eig1c), kind=8) !original wrong here
            eig3r = aaaa + bbbb - aa / 3.0d0

            rortexEigr(i,j,k) = abs(eig3r)

            ! real right eigenvector

            delta1 = (a(1,1)-eig3r)*(a(2,2)-eig3r) - a(2,1)*a(1,2)
            delta2 = (a(2,2)-eig3r)*(a(3,3)-eig3r) - a(2,3)*a(3,2)
            delta3 = (a(1,1)-eig3r)*(a(3,3)-eig3r) - a(1,3)*a(3,1)

            if(delta1 == 0.0d0 .and. delta2 == 0.0d0 .and. delta3 == 0.0d0) then
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
              vr(3) = 1.0d0

            else if(abs(delta2) >= abs(delta1) .and.                           &
                    abs(delta2) >= abs(delta3)) then

              vr(1) = 1.0d0
              vr(2) = (-(a(3,3)-eig3r)*a(2,1) +         a(2,3)*a(3,1))/delta2
              vr(3) = (         a(3,2)*a(2,1) - (a(2,2)-eig3r)*a(3,1))/delta2

            else if(abs(delta3) >= abs(delta1) .and.                           &
                    abs(delta3) >= abs(delta2)) then

               vr(1) = (-(a(3,3)-eig3r)*a(1,2) +         a(1,3)*a(3,2))/delta3
               vr(2) = 1.0d0
               vr(3) = (         a(3,1)*a(1,2) - (a(1,1)-eig3r)*a(3,2))/delta3

            else

              write(*,*) 'ERROR: '
              write(*,*) delta1, delta2, delta3
              stop

            end if

            temp = sqrt(vr(1)**2 + vr(2)**2 + vr(3)**2)

            vr(1) = vr(1) / temp
            vr(2) = vr(2) / temp
            vr(3) = vr(3) / temp

            call rotation(z0, vr, qqq)

            vg = matmul(transpose(qqq), a)
            vg = matmul(vg, qqq)

            alpha = 0.5*sqrt((vg(2,2)-vg(1,1))**2+(vg(2,1)+vg(1,2))**2)
            beta  = 0.5*(vg(2,1)-vg(1,2))

            if(beta**2 > alpha**2) then

              if(beta > 0.0) then
                rm = 2*(beta-alpha)
                rortex_x(i, j, k) = rm*vr(1)
                rortex_y(i, j, k) = rm*vr(2)
                rortex_z(i, j, k) = rm*vr(3)
              else
                rm = 2*(beta+alpha)
                rortex_x(i, j, k) = rm*vr(1)
                rortex_y(i, j, k) = rm*vr(2)
                rortex_z(i, j, k) = rm*vr(3)
              end if

            else

              rortex_x(i, j, k) = 0.0
              rortex_y(i, j, k) = 0.0
              rortex_z(i, j, k) = 0.0

            end if

            rortex_mag(i,j,k) = sqrt(rortex_x(i,j,k)**2                        &
                                    +rortex_y(i,j,k)**2                        &
                                    +rortex_z(i,j,k)**2)

            lamb_ci(i,j,k) = abs(aimag(eig1c))

          else ! three real roots

            rortex_x(i,j,k) = 0.0d0
            rortex_y(i,j,k) = 0.0d0
            rortex_z(i,j,k) = 0.0d0
            rortex_mag(i,j,k) = 0.0d0

          end if

          ljmtmp2=sqrt(rortex_mag_x(i,j,k)**2 + rortex_mag_y(i,j,k)**2 &
                +rortex_mag_z(i,j,k)**2)

          ljmtmp=sqrt((rortex_mag_y(i,j,k)*rortex_z(i,j,k)-rortex_mag_z(i,j,k)*rortex_y(i,j,k))**2 &
              +(rortex_mag_z(i,j,k)*rortex_x(i,j,k)-rortex_mag_x(i,j,k)*rortex_z(i,j,k))**2 &
              +(rortex_mag_x(i,j,k)*rortex_y(i,j,k)-rortex_mag_y(i,j,k)*rortex_x(i,j,k))**2)

          if (ljmtmp2>1.0e-6) then
              rortexEigr(i,j,k)=abs(rortex_mag_x(i,j,k)*vr(1)+rortex_mag_y(i,j,k)*vr(2)+rortex_mag_z(i,j,k)*vr(3) )/ljmtmp2
          else
              rortexEigr(i,j,k)=0.0
          end if

         omgGradXRortex(i,j,k)=1.0

         if(localOmega(i,j,k) < 0.51) then
             omgGradXRortex(i,j,k)=0.0
         end if

         rortex_mag_x(i,j,k) = rortex_mag_x(i,j,k)*omgGradXRortex(i,j,k)
         rortex_mag_y(i,j,k) = rortex_mag_y(i,j,k)*omgGradXRortex(i,j,k)
         rortex_mag_z(i,j,k) = rortex_mag_z(i,j,k)*omgGradXRortex(i,j,k)

        end do
      end do
    end do

    call cpu_time(t6)

    write(*,*) 'calculation time: ', t6 - t5

    1234 continue


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

    do k = 1, kmax_local
      do j = 1, jmax_local
        do i = 1, imax_local

          a(1,1) = dudx(i, j, k)
          a(1,2) = dudy(i, j, k)
          a(1,3) = dudz(i, j, k)
          a(2,1) = dvdx(i, j, k)
          a(2,2) = dvdy(i, j, k)
          a(2,3) = dvdz(i, j, k)
          a(3,1) = dwdx(i, j, k)
          a(3,2) = dwdy(i, j, k)
          a(3,3) = dwdz(i, j, k)

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
          delta = -delta / 108.d0

          if(delta > 0.0d0) then ! one real root and two complex conjugate roots

            aaaa = -sign(1.0, rr)*(abs(rr)+sqrt(delta))**(1.0d0/3.0d0)

            if(aaaa == 0.0d0) then
              bbbb = 0.0d0
            else
              bbbb = qq / aaaa
            end if

            eig1c = cmplx(-0.5d0*(aaaa+bbbb)-aa/3.0, 0.5*sqrt(3.0)*(aaaa-bbbb), kind=8)
            eig2c = cmplx(real(eig1c), -aimag(eig1c), kind=8)
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

    write(*,*) 'Finished Modified Liutex Omega'

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


    outputfilename1 = 'data/'//trim(outputfilename)//'_'//chars//".xyz"
    outputfilename2 = 'data/'//trim(outputfilename)//'_'//chars//".fun"

    ! write coordinates
    open(fin1, file=trim(outputfilename1), form='unformatted', action='write')
    write(fin1) imax_local, jmax_local, kmax_local
    write(fin1)                                                                &
               (((x(i, j, k), i=istart,iend), j=jstart,jend), k=kstart,kend),  &
               (((y(i, j, k), i=istart,iend), j=jstart,jend), k=kstart,kend),  &
               (((z(i, j, k), i=istart,iend), j=jstart,jend), k=kstart,kend)

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
    f(:,:,:,7) = localOmega
    f(:,:,:,8) = rortex_x
    f(:,:,:,9) = rortex_y
    f(:,:,:,10) = rortex_z
    f(:,:,:,11) = rortex_mag
    f(:,:,:,12) = Qmethod
    f(:,:,:,13) = isLocExtrOmg
    f(:,:,:,14) = rortex_mag_x
    f(:,:,:,15) = rortex_mag_y
    f(:,:,:,16) = rortex_mag_z
    f(:,:,:,17) = omega_l

    open(fin1, file=trim(outputfilename2), form='unformatted', action='write')

    write(fin1) imax_local, jmax_local, kmax_local, nvar

    write(fin1) ((((f(i, j, k, nx), i=1,imax_local), &
                                    j=1,jmax_local), &
                                    k=1,kmax_local), &
                                    nx=1,nvar)

    close(fin1)

    deallocate(f)

    write(*,*) trim(outputfilename), ' finish successfully'

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

    deallocate(omega)

    deallocate(u)
    deallocate(v)
    deallocate(w)

    deallocate(rortex_x)
    deallocate(rortex_y)
    deallocate(rortex_z)

    deallocate(rortex_mag)
    deallocate(Qmethod)
    deallocate(lamb_ci)
    deallocate(localOmega)

    deallocate(QGlob)
    deallocate(QLoc_x)
    deallocate(QLoc_y)
    deallocate(QLoc_z)
    deallocate(omgGradXRortex)
    deallocate(omgLoc_x)
    deallocate(omgLoc_y)
    deallocate(omgLoc_z)
    deallocate(isLocExtrOmg)
    deallocate(rortexEigr)
    deallocate(rortex_mag_x)
    deallocate(rortex_mag_y)
    deallocate(rortex_mag_z)
    deallocate(rortexGlob_mag)

    deallocate(o_alpha)
    deallocate(o_beta)

    deallocate(lambda_cr)
    deallocate(lambda_r)

    deallocate(omega_l)

  end do

end program DNSUTA_extract



! calculate the rortex
! a velocity gradient tensor, vor: vorticity,vr: rotational axis, rorMag:rotational strength
subroutine cal_rortex(a, vor, vr,rorMag)
  implicit none

  real(8), intent(in) :: a(3,3)
  real(8), intent(in) :: vor(3)
  real(8), intent(out) :: vr(3), rorMag

  real(8):: aa,bb,cc,delta,rr,aaaa,bbbb,qq,delta1,delta2,delta3,temp
  complex(8):: eig1c,eig2c,eig3r
  real(8) :: tt(3,3)
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
          vr=0.0
          rorMag=0.0

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

end subroutine cal_rortex


subroutine rotation(u, v, r)
!-------------------------------------------------------------------------------
! calculate rotation matrix r which rotates unit vector u to unit vector v
!-------------------------------------------------------------------------------

  implicit none

  real(8), intent(in) :: u(3)
  real(8), intent(in) :: v(3)
  real(8), intent(out) :: r(3, 3)

  real(8) :: a(3)
  real(8) :: aa
  real(8) :: t
  real(8) :: alpha
  real(8) :: c, s

  real(8), parameter :: eps = 1.0e-10

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

