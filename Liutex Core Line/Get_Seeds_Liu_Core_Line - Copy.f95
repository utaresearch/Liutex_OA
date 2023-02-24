program get_seeds_liu_core_line
!-------------------------------------------------------------------------------
! This program reads DNSUTA_extractOA.f95 grid (.xyz) and
! function (.fun) files and finds the seed points for the Liutex Core Line.
!
! Author: Oscar Alvarez (oscar.alvarez@mavs.uta.edu)
! Department of Mathematics, University of Texas at Arlington,
! Arlington, Texas, USA
! Spring 2023
!-------------------------------------------------------------------------------

  implicit none

  real(8), parameter :: pi = 4.0d0*atan(1.0d0)
  
  real(8), dimension(:,:,:,:), allocatable :: f
  real(8), dimension(:,:,:,:), allocatable :: liu_mag_grad
  
  real(8), dimension(:,:,:), allocatable :: liutex_x, liutex_y, liutex_z, liutex_mag
  real(8), dimension(:,:,:), allocatable :: xi_x, xi_y, xi_z
  real(8), dimension(:,:,:), allocatable :: eta_x, eta_y, eta_z
  real(8), dimension(:,:,:), allocatable :: zeta_x, zeta_y, zeta_z
  real(8), dimension(:,:,:), allocatable :: dmol_dx, dmol_dy, dmol_dz
  real(8), dimension(:,:,:), allocatable :: mol_x, mol_y, mol_z
  real(8), dimension(:,:,:), allocatable :: x, y, z
  real(8), dimension(:,:,:), allocatable :: u, v, w
  real(8), dimension(:,:,:), allocatable :: mol

  real(8) :: x_xi, x_eta, x_zeta
  real(8) :: y_xi, y_eta, y_zeta
  real(8) :: z_xi, z_eta, z_zeta
  real(8) :: det
  real(8) :: mol_xi, mol_eta, mol_zeta
  real(8) :: cross_mag
  real(8) :: cross1, cross2, cross3
 
  integer, dimension(:,:,:), allocatable :: is_seed
  integer, dimension(:,:,:), allocatable :: loc_max
  
  integer :: s_counter
  integer :: i, j, k, nx
  integer :: imax, jmax, kmax
  
  !! file handling
  integer, parameter :: fin1 = 10, fin2 = 20, fout1 = 30
  character(100) :: inputfilename
  character(100) :: gridfilename, funcfilename
  character(100) :: seedsfilename, outputfilename
  character(100) :: msg
  character(100) :: ignore, datafileprefix
  character(6) :: chars
  integer :: ios
  integer :: nvar
  logical :: lexist
  integer :: f_start, f_end
  integer :: skip
  integer :: istart, iend, jstart, jend, kstart, kend
  integer :: iii
  integer :: imax_local, jmax_local, kmax_local


  inputfilename = 'input.txt'

  !! inquire whether the input file exists
  inquire(file=inputfilename, exist=lexist)
  if(.not. lexist) then
    write(*,*) 'ERROR: no file '//inputfilename
    stop
  end if

  !! read input file
  open(fin1, file=inputfilename, form='formatted', action='read')
  read(fin1, *) f_start
  read(fin1, *) f_end
  read(fin1, *) skip
  read(fin1, *) istart, iend
  read(fin1, *) jstart, jend
  read(fin1, *) kstart, kend
  read(fin1, *) ignore
  read(fin1, *) datafileprefix

  close(fin1)


  imax_local = iend - istart + 1
  jmax_local = jend - jstart + 1
  kmax_local = kend - kstart + 1

  call itoa6(f_start, chars)

  gridfilename = trim(datafileprefix)//'_'//chars//'.xyz'

  !! check if the grid file exists
  inquire(file=trim(gridfilename), exist=lexist)
  if(.not. lexist) then
    write(*,*) 'ERROR: no grid file ', trim(gridfilename)
    stop
  end if

  write(*,*)
  write(*,*) 'reading ', trim(gridfilename)

  !! open grid file
  open(fin1, file=trim(gridfilename), form='unformatted', action='read', iostat=ios, iomsg=msg)
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

  !! read coordinates
  read(fin1, iostat=ios, iomsg=msg) (((x(i, j, k), i=1,imax), j=1,jmax), k=1,kmax),  &
                                    (((y(i, j, k), i=1,imax), j=1,jmax), k=1,kmax),  &
                                    (((z(i, j, k), i=1,imax), j=1,jmax), k=1,kmax)

  close(fin1)

  print*, 'Grid File Read Complete !!'

  do iii = f_start, f_end, skip

    call itoa6(iii, chars)

    funcfilename = trim(datafileprefix)//'_'//chars//'.fun'

    !! check if the function file exists
    inquire(file=trim(funcfilename), exist=lexist)
    if(.not. lexist) then
      write(*,*) 'ERROR: no data file ', funcfilename
      stop
    end if

    write(*,*)
    write(*,*) 'reading ', trim(funcfilename)

    !! open function file
    open(fin2, file=trim(funcfilename), form='unformatted', action='read', iostat=ios, iomsg=msg)
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

    ! Reading Data
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

    read(fin2, iostat=ios, iomsg=msg) ((((f(i, j, k, nx), i=1,imax), j=1,jmax), k=1,kmax), nx=1,nvar)
    if(ios /= 0) then
      write(*,*) 'ERROR: ', msg
      stop
    end if

    close(fin2)

    print*, "Function File Read Complete !!"

    allocate(u(imax, jmax, kmax))
    allocate(v(imax, jmax, kmax))
    allocate(w(imax, jmax, kmax))

    allocate(liutex_x(imax,jmax,kmax))
    allocate(liutex_y(imax,jmax,kmax))
    allocate(liutex_z(imax,jmax,kmax))

    allocate(liutex_mag(imax,jmax,kmax))
    allocate(liu_mag_grad(imax,jmax,kmax,3))

    allocate(mol(imax, jmax, kmax))

    u = f(:,:,:,1)
    v = f(:,:,:,2)
    w = f(:,:,:,3)

    liutex_x = f(:,:,:,8)
    liutex_y = f(:,:,:,9)
    liutex_z = f(:,:,:,10)
    liutex_mag = f(:,:,:,11)

    liu_mag_grad(:,:,:,1) = f(:,:,:,14) ! gradient in x direction
    liu_mag_grad(:,:,:,2) = f(:,:,:,15) ! gradient in y direction
    liu_mag_grad(:,:,:,3) = f(:,:,:,16) ! gradient in z direction

    mol = f(:,:,:,17)

    deallocate(f)

    write(*,*) 'Calculating Liutex Core Line Seed Points'

    allocate(xi_x(imax, jmax, kmax))
    allocate(xi_y(imax, jmax, kmax))
    allocate(xi_z(imax, jmax, kmax))

    allocate(eta_x(imax, jmax, kmax))
    allocate(eta_y(imax, jmax, kmax))
    allocate(eta_z(imax, jmax, kmax))

    allocate(zeta_x(imax, jmax, kmax))
    allocate(zeta_y(imax, jmax, kmax))
    allocate(zeta_z(imax, jmax, kmax))

    allocate(dmol_dx(imax, jmax, kmax))
    allocate(dmol_dy(imax, jmax, kmax))
    allocate(dmol_dz(imax, jmax, kmax))

    !! Calculating gradient of Modified omega liutex (mol)
    do k = 1, kmax
      do j = 1, jmax
        do i = 1, imax

          !! First condition of finding core line
          if(mol(i,j,k) >= 0.52d0) then

            !! Finite-Difference for first partial derivatives (gradient)
            if(i == 1) then
              mol_xi = mol(2, j, k) - mol(1, j, k)

              x_xi = x(2, j, k) - x(1, j, k)
              y_xi = y(2, j, k) - y(1, j, k)
              z_xi = z(2, j, k) - z(1, j, k)
            else if(i == imax) then
              mol_xi = mol(imax, j, k) - mol(imax-1, j, k)

              x_xi = x(imax, j, k)-x(imax-1, j, k)
              y_xi = y(imax, j, k)-y(imax-1, j, k)
              z_xi = z(imax, j, k)-z(imax-1, j, k)
            else
              mol_xi = 0.5d0*(mol(i+1, j, k) - mol(i-1, j, k))

              x_xi = 0.5d0*(x(i+1, j, k) - x(i-1, j, k))
              y_xi = 0.5d0*(y(i+1, j, k) - y(i-1, j, k))
              z_xi = 0.5d0*(z(i+1, j, k) - z(i-1, j, k))
            end if

            if(j == 1) then
              mol_eta = mol(i, 2, k) - mol(i, 1, k)

              x_eta = x(i, 2, k) - x(i, 1, k)
              y_eta = y(i, 2, k) - y(i, 1, k)
              z_eta = z(i, 2, k) - z(i, 1, k)
            else if(j == jmax) then
              mol_eta = mol(i, jmax, k) - mol(i, jmax-1, k)

              x_eta = x(i, jmax, k) - x(i, jmax-1, k)
              y_eta = y(i, jmax, k) - y(i, jmax-1, k)
              z_eta = z(i, jmax, k) - z(i, jmax-1, k)
            else
              mol_eta = 0.5d0*(mol(i, j+1, k) - mol(i, j-1, k))

              x_eta = 0.5d0*(x(i, j+1, k) - x(i, j-1, k))
              y_eta = 0.5d0*(y(i, j+1, k) - y(i, j-1, k))
              z_eta = 0.5d0*(z(i, j+1, k) - z(i, j-1, k))
            end if

            if(k == 1) then
              mol_zeta = mol(i, j, 2) - mol(i, j, 1)
              x_zeta = x(i, j, 2) - x(i, j, 1)
              y_zeta = y(i, j, 2) - y(i, j, 1)
              z_zeta = z(i, j, 2) - z(i, j, 1)
            else if(k == kmax) then
              mol_zeta = mol(i, j, kmax) - mol(i, j, kmax-1)

              x_zeta = x(i, j, kmax) - x(i, j, kmax-1)
              y_zeta = y(i, j, kmax) - y(i, j, kmax-1)
              z_zeta = z(i, j, kmax) - z(i, j, kmax-1)
            else
              mol_zeta = 0.5d0*(mol(i, j, k+1) - mol(i, j, k-1))

              x_zeta = 0.5d0*(x(i, j, k+1) - x(i, j, k-1))
              y_zeta = 0.5d0*(y(i, j, k+1) - y(i, j, k-1))
              z_zeta = 0.5d0*(z(i, j, k+1) - z(i, j, k-1))
            end if

            !! determinant of Jacobian
            det =   x_xi*(y_eta*z_zeta-y_zeta*z_eta)  &
                  - x_eta*(y_xi*z_zeta-y_zeta*z_xi)   &
                  + x_zeta*(y_xi*z_eta-y_eta*z_xi)

            det = 1.0d0 / det

            xi_x(i,j,k) = det*(y_eta*z_zeta - y_zeta*z_eta)
            xi_y(i,j,k) = det*(x_zeta*z_eta - x_eta*z_zeta)
            xi_z(i,j,k) = det*(x_eta*y_zeta - x_zeta*y_eta)

            eta_x(i,j,k) = det*(y_zeta*z_xi - y_xi*z_zeta)
            eta_y(i,j,k) = det*(x_xi*z_zeta - x_zeta*z_xi)
            eta_z(i,j,k) = det*(x_zeta*y_xi - x_xi*y_zeta)

            zeta_x(i,j,k) = det*(y_xi*z_eta - y_eta*z_xi)
            zeta_y(i,j,k) = det*(x_eta*z_xi - x_xi*z_eta)
            zeta_z(i,j,k) = det*(x_xi*y_eta - x_eta*y_xi)

            !! First partial derivatives (Gradient of Modified-Omega-Liutex)
            dmol_dx(i,j,k) = mol_xi*xi_x(i,j,k) + mol_eta*eta_x(i,j,k) + mol_zeta*zeta_x(i,j,k)
            dmol_dy(i,j,k) = mol_xi*xi_y(i,j,k) + mol_eta*eta_y(i,j,k) + mol_zeta*zeta_y(i,j,k)
            dmol_dz(i,j,k) = mol_xi*xi_z(i,j,k) + mol_eta*eta_z(i,j,k) + mol_zeta*zeta_z(i,j,k)

          end if

        end do
      end do
    end do
   
    allocate(is_seed(imax,jmax,kmax))

    seedsfilename = 'seeds_'//trim(datafileprefix)//'_'//chars//'.txt'
    open(fout1, file=trim(seedsfilename), form='formatted', action='write')

    s_counter = 0   !! Seed point counter

    do k = 1, kmax
      do j = 1, jmax
        do i = 1, imax

          if(mol(i,j,k) >= 0.52d0) then 
            !-------------------------------------------------------------------
            !            Finding (gradient of Liutex) x r = 0
            !
            ! This is the definition of the Liutex Core Line.
            !-------------------------------------------------------------------

            !  cross1 = liu_mag_grad(i,j,k,2)*liutex_z(i,j,k)                 &
            !              - liu_mag_grad(i,j,k,3)*liutex_y(i,j,k)
            !  cross2 = liu_mag_grad(i,j,k,3)*liutex_x(i,j,k)                 &
            !              - liu_mag_grad(i,j,k,1)*liutex_z(i,j,k)
            !  cross3 = liu_mag_grad(i,j,k,1)*liutex_y(i,j,k)                 &
            !              - liu_mag_grad(i,j,k,2)*liutex_x(i,j,k)

            !-------------------------------------------------------------------
            !            Finding (mol_vec) x r = 0
            !-------------------------------------------------------------------

            cross1 = mol_y(i,j,k)*liutex_z(i,j,k) - mol_z(i,j,k)*liutex_y(i,j,k)
            cross2 = mol_z(i,j,k)*liutex_x(i,j,k) - mol_x(i,j,k)*liutex_z(i,j,k)
            cross3 = mol_x(i,j,k)*liutex_y(i,j,k) - mol_y(i,j,k)*liutex_x(i,j,k)

            cross_mag = sqrt(cross1*cross1 + cross2*cross2 + cross3*cross3)

            if (cross_mag == 0.0) then
              is_seed(i,j,k) = 1
              s_counter = s_counter + 1

              !! write x, y, z seed point location to file
              write(fout1,*) x(i,j,k), y(i,j,k), z(i,j,k)
            end if
          
          end if

        end do
      end do
    end do

    close(fout1)

    deallocate(dmol_dx)
    deallocate(dmol_dy)
    deallocate(dmol_dz)
    deallocate(loc_max)

    deallocate(u)
    deallocate(v)
    deallocate(w)

    deallocate(liutex_x)
    deallocate(liutex_y)
    deallocate(liutex_z)

    deallocate(liutex_mag)
    deallocate(liu_mag_grad)

    deallocate(mol)

    !! Check if any seed points were found
    write(*,*) ''

    if(s_counter == 0) then
      write(*,*) '!! NO SEED POINTS FOUND: '//trim(datafileprefix)//'_'//chars//' !!'
    else
      write(*,*) 'Finished. Located ', s_counter, ' seed points: '//trim(datafileprefix)//'_'//chars
    end if

    write(*,*) ''

  end do

  deallocate(x)
  deallocate(y)
  deallocate(z)

  deallocate(xi_x)
  deallocate(xi_y)
  deallocate(xi_z)

  deallocate(eta_x)
  deallocate(eta_y)
  deallocate(eta_z)

  deallocate(zeta_x)
  deallocate(zeta_y)
  deallocate(zeta_z)

end program get_seeds_liu_core_line


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
