program main
!!!=================================================================================================
!!! Liutex Core Line generation program.
!!! 
!!! First compile module:  then compile together, i.e.:
!!! > gfortran -c liutex_mod.f03
!!! 
!!! Then compile this program:
!!! > gfortran -c coreline.f03
!!!
!!! Then compile the output files (.o) together:
!!! > gfortran liutex_mod.o coreline.o -o coreline.exe
!!! 
!!! By: Oscar Alvarez
!!!=================================================================================================
    use liutex_mod

    implicit none
    
    !! Parameter values
    real(8), parameter :: pi = 4.d0*atan(1.d0)

    !! File handling
    integer, parameter :: fin1 = 10, fin2 = 20, fout1 = 30
    character(100) :: inputfilename, gridfilename, funcfilename, outputfilename, pointfilename
    character(100) :: msg, ignore, datafileprefix
    character(6) :: chars
    integer :: ios, nvar
    integer :: f_start, f_end, skip
    integer :: istart, iend, jstart, jend, kstart, kend
    integer :: iii
    logical :: lexist

    !! Other stuff
    real(8), dimension(:,:,:,:), allocatable :: f, f2
    real(8), dimension(:,:,:,:), allocatable :: r, lmg
    real(8), dimension(:,:,:), allocatable :: x, y, z
    real(8), dimension(:,:,:), allocatable :: r_mag, omega_l
    real(8), dimension(:,:,:), allocatable :: liutex_core_x, liutex_core_y, liutex_core_z
    real(8), dimension(:,:,:), allocatable :: div_lmg, div_r
    real(8), dimension(3) :: lmg_vec, r_vec
    real(8) :: lmg_norm, r_norm
    real(8) :: omega_l_tol, dot_tol, div_r_tol, div_lmg_tol
    logical :: omega_l_condition, div_lmg_condition, div_r_condition, conditions
    integer :: imax, jmax, kmax, nx
    integer :: i, j, k
    integer :: liutex_core_counter

    
    !! Program start
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

    write(chars,"(I6)") f_start
    chars = trim(chars)
    
    gridfilename = 'data/'//trim(datafileprefix)//'_'//chars//'.xyz'

    !! Check if the grid file exists
    inquire(file=trim(gridfilename), exist=lexist)
    if(.not. lexist) then
        write(*,*) 'ERROR: no grid file ', trim(gridfilename)
        stop
    end if

    write(*,*) 'reading ', trim(gridfilename)

    !! Open grid file
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
    write(*,*) 'Grid file (.xyz) dimensions:'
    write(*,*) 'imax = ', imax, ' jmax = ', jmax, ' kmax = ', kmax
    write(*,*)
    !! read coordinates
    read(fin1, iostat=ios, iomsg=msg) (((x(i, j, k), i=1,imax), j=1,jmax), k=1,kmax),  &
                                      (((y(i, j, k), i=1,imax), j=1,jmax), k=1,kmax),  &
                                      (((z(i, j, k), i=1,imax), j=1,jmax), k=1,kmax)

    close(fin1)

    write(*,*) 'Grid File Read Complete !!'

    do iii = f_start, f_end, skip

        write(chars, "(I6)") iii

        funcfilename = 'data/'//trim(datafileprefix)//'_'//chars//'.fun'

        !! Check if the function file exists
        inquire(file=trim(funcfilename), exist=lexist)
        if(.not. lexist) then
            write(*,*) 'ERROR: no data file ', funcfilename
            stop
        end if

        write(*,*) 'reading ', funcfilename

        !! Open function file
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

        allocate(r(imax,jmax,kmax,3))
        allocate(r_mag(imax,jmax,kmax))
        allocate(lmg(imax,jmax,kmax,3))
        allocate(omega_l(imax,jmax,kmax))

        r      = f(:,:,:, 8:10)
        r_mag  = f(:,:,:, 11)
        lmg    = f(:,:,:, 14:16)
        omega_l= f(:,:,:, 17)
        
        allocate(liutex_core_x(imax,jmax,kmax))   !! New Liutex Core Vector
        allocate(liutex_core_y(imax,jmax,kmax))   !! New Liutex Core Vector
        allocate(liutex_core_z(imax,jmax,kmax))   !! New Liutex Core Vector
        
        !! Normalize vector fields
        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax
                    !! Create vectors.
                    lmg_vec = (/ lmg(i,j,k,1), lmg(i,j,k,2), lmg(i,j,k,3) /)
                    r_vec = (/ r(i,j,k,1), r(i,j,k,2), r(i,j,k,3) /)

                    lmg_norm = norm2(lmg_vec)
                    r_norm = norm2(r_vec)

                    !! Check if the norm is zero and if not then divide.
                    if (lmg_norm .ne. 0.d0) then
                        lmg(i,j,k,1) = lmg(i,j,k,1) / lmg_norm
                        lmg(i,j,k,2) = lmg(i,j,k,2) / lmg_norm
                        lmg(i,j,k,3) = lmg(i,j,k,3) / lmg_norm
                    end if

                    if (r_norm .ne. 0.d0) then
                        r(i,j,k,1) = r(i,j,k,1) / r_norm
                        r(i,j,k,2) = r(i,j,k,2) / r_norm
                        r(i,j,k,3) = r(i,j,k,3) / r_norm
                    end if 

                end do
            end do 
        end do
        
        write(*,*) 'Finding the Liutex Core Vector Field'

        pointfilename = 'data/corepoints_'//trim(datafileprefix)//'_'//chars//'.dat'
        open(fout1, file=trim(pointfilename), form='formatted', action='write')

        liutex_core_x = 0.d0
        liutex_core_y = 0.d0
        liutex_core_z = 0.d0

        allocate(div_lmg(imax,jmax,kmax))
        allocate(div_r(imax,jmax,kmax))

        !! Finding the divergence of vectors.
        div_lmg = divergence(abs(lmg(:,:,:,1)), abs(lmg(:,:,:,2)), abs(lmg(:,:,:,3)), x, y, z, imax, jmax, kmax)
        ! div_r   = divergence(abs(r(:,:,:,1)), abs(r(:,:,:,2)), abs(r(:,:,:,3)), x, y, z, imax, jmax, kmax)

        liutex_core_counter = 0

        !! Apply conditions for liutex core vector

        omega_l_tol = 0.51d0        !! tolerance for omega liutex
        dot_tol     = 1.d-6         !! tolerance for dot product
        div_r_tol   = 1.d-6         !! tolerance for divergence of liutex condition
        div_lmg_tol = 1.d-4         !! tolerance for divergence of liutex gradient condition

        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax

                    omega_l_condition = (omega_l(i,j,k) > omega_l_tol)

                    if (omega_l_condition) then

                        ! div_r_condition = (abs(div_r(i,j,k)) < div_r_tol)
                        div_lmg_condition = (abs(div_lmg(i,j,k)) > div_lmg_tol)

                        ! conditions = div_r_condition
                        conditions = div_lmg_condition
                        ! conditions = div_lmg_condition .and. div_r_condition

                        if (conditions) then

                            !! Liutex Core Line (liutex_core) vector components
                            liutex_core_x(i,j,k) = r_mag(i,j,k) * lmg(i,j,k,1)
                            liutex_core_y(i,j,k) = r_mag(i,j,k) * lmg(i,j,k,2)
                            liutex_core_z(i,j,k) = r_mag(i,j,k) * lmg(i,j,k,3)

                            liutex_core_counter = liutex_core_counter + 1

                            !! Write liutex core points to file
                            ! write(fout1,*) x(i,j,k), y(i,j,k), z(i,j,k)

                        end if

                    end if
        
                end do
            end do
        end do

        close(fout1)

        !! Check if any liutex core lines were found.
        if (liutex_core_counter == 0) then
            !! This should ONLY happen if there there is no turbulence/vortices.
            write(*,*)
            write(*,*) 'NO LIUTEX CORE LINES FOUND FOR: ', funcfilename
            goto 100
        else 
            write(*,*)
            write(*,*) 'Liutex Core Lines Found: ', liutex_core_counter
            write(*,*)
        end if

        write(*,*) 'Writing Function (.fun) file.'

        nvar = nvar + 3
        allocate(f2(imax,jmax,kmax,nvar))

        f2(:,:,:, 1:nvar-3) = f
        deallocate(f)
        
        f2(:,:,:, nvar-2)   = liutex_core_x
        f2(:,:,:, nvar-1)   = liutex_core_y
        f2(:,:,:, nvar)     = liutex_core_z
        
        outputfilename = 'data/corelines_'//trim(datafileprefix)//'_'//chars//'.fun'

        open(fout1, file=trim(outputfilename), form='unformatted', action='write')
        write(fout1) imax, jmax, kmax, nvar
        write(fout1) ((((f2(i, j, k, nx), i=1,imax), j=1,jmax), k=1,kmax), nx=1,nvar)
        close(fout1)
        
        write(*,*) 'New function (.fun) file created: ', outputfilename

        100 continue

        deallocate(r)
        deallocate(r_mag)
        deallocate(lmg)
        deallocate(omega_l)
        deallocate(liutex_core_x)
        deallocate(liutex_core_y)
        deallocate(liutex_core_z)
        
    end do

    write(*,*) 'Program complete.'
  
end program main