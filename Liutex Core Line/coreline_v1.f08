program main
    !!!=================================================================================================
    !!! Liutex Core Line generation program.
    !!! 
    !!! First compile module then compile final program, e.g.:
    !!! > gfortran -c liutex_mod.f03
    !!! > gfortran coreline.f03 -o coreline.exe
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
        real(8), dimension(:,:,:,:), allocatable :: r, l_mag_gradient
        real(8), dimension(:,:,:), allocatable :: x, y, z
        real(8), dimension(:,:,:), allocatable :: l_mag, omega_l
        real(8), dimension(:,:,:), allocatable :: l_core_x, l_core_y, l_core_z
        real(8), dimension(3) :: lmg_vec, r_vec
        real(8) :: lmg_norm, r_norm
        real(8) :: tol1, tol2, tol3
        logical :: condition1, condition2, condition3, all_conditions
        integer :: imax, jmax, kmax, nx
        integer :: i, j, k
        integer :: l_core_counter
    
        
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
    
            write(chars, "(I6)") iii
    
            funcfilename = 'data/'//trim(datafileprefix)//'_'//chars//'.fun'
    
            !! Check if the function file exists
            inquire(file=trim(funcfilename), exist=lexist)
            if(.not. lexist) then
                write(*,*) 'ERROR: no data file ', funcfilename
                stop
            end if
    
            write(*,*)
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
            allocate(l_mag(imax,jmax,kmax))
            allocate(l_mag_gradient(imax,jmax,kmax,3))
            allocate(omega_l(imax,jmax,kmax))
    
            r(:,:,:,:) = f(:,:,:, 8:10)
            l_mag = f(:,:,:, 11)
            l_mag_gradient(:,:,:,:) = f(:,:,:, 14:16)
            omega_l = f(:,:,:, 17)
            
            allocate(l_core_x(imax,jmax,kmax))   !! New Liutex Core Vector
            allocate(l_core_y(imax,jmax,kmax))   !! New Liutex Core Vector
            allocate(l_core_z(imax,jmax,kmax))   !! New Liutex Core Vector
            
            write(*,*) 'Finding the Liutex Core Vector Field'
            
            pointfilename = 'data/corepoints_'//trim(datafileprefix)//'_'//chars//'.dat'
            open(fout1, file=trim(pointfilename), form='formatted', action='write')
    
            l_core_x = 0.d0
            l_core_y = 0.d0
            l_core_z = 0.d0
    
            !! Apply conditions for liutex core vector
    
            tol1 = 0.51d0   !! tolerance for condition 1
            tol2 = 1.d-6    !! tolerance for condition 2
            tol3 = 1.d-3    !! tolerance for condition 3
    
            l_core_counter = 0
    
            do k = 1, kmax
                do j = 1, jmax
                    do i = 1, imax
    
                        condition1 = (omega_l(i,j,k) >= tol1)
                        
                        if (condition1) then
    
                            lmg_vec = (/ l_mag_gradient(i,j,k,1), l_mag_gradient(i,j,k,2), l_mag_gradient(i,j,k,3) /)
                            r_vec = (/ r(i,j,k,1), r(i,j,k,2), r(i,j,k,3) /)
    
                            lmg_norm = norm2(lmg_vec)
                            r_norm = norm2(r_vec)
    
                            if (lmg_norm .ne. 0.d0) then
                                lmg_vec = lmg_vec / lmg_norm
                            end if
    
                            if (r_norm .ne. 0.d0) then
                                r_vec = r_vec / r_norm
                            end if
    
                            !! Conditions
                            condition2 = (lmg_norm <= tol2)
                            
                            condition3 = (norm2(cross_product_3d(r_vec, lmg_vec)) <= tol3)
                            
                            all_conditions = condition2 .and. condition3
    
                            if (all_conditions) then
    
                                !! Liutex Core Line (l_core) vector components
                                l_core_x(i,j,k) = l_mag(i,j,k) * lmg_vec(1)
                                l_core_y(i,j,k) = l_mag(i,j,k) * lmg_vec(2)
                                l_core_z(i,j,k) = l_mag(i,j,k) * lmg_vec(3)
    
                                l_core_counter = l_core_counter + 1
    
                                !! Write liutex core points to file
                                write(fout1,*) x(i,j,k), y(i,j,k), z(i,j,k)
    
                            end if
    
                        end if
            
                    end do
                end do
            end do
    
            close(fout1)
    
            !! Check if any liutex core lines were found.
            if (l_core_counter == 0) then
                !! This should ONLY happen if there there is no turbulence/vortices.
                write(*,*)
                write(*,*) 'NO LIUTEX CORE LINES FOUND FOR: ', funcfilename
                goto 1
            else 
                write(*,*)
                write(*,*) 'Liutex Core Lines Found: ', l_core_counter
            end if
    
            nvar = nvar + 3
            allocate(f2(imax,jmax,kmax,nvar))
    
            f2(:,:,:, 1:nvar-3) = f
            deallocate(f)
            
            f2(:,:,:, nvar-2)   = l_core_x
            f2(:,:,:, nvar-1)   = l_core_y
            f2(:,:,:, nvar)     = l_core_z
            
            write(*,*) 'Writing Function (.fun) file.'
    
            outputfilename = 'data/corelines_'//trim(datafileprefix)//'_'//chars//'.fun'
    
            open(fout1, file=trim(outputfilename), form='unformatted', action='write')
            write(fout1) imax, jmax, kmax, nvar
            write(fout1) ((((f2(i, j, k, nx), i=1,imax), j=1,jmax), k=1,kmax), nx=1,nvar)
            close(fout1)
            
            write(*,*)
            write(*,*) 'New function (.fun) file created: ', outputfilename
    
            1 continue
    
            deallocate(r)
            deallocate(l_mag)
            deallocate(l_mag_gradient)
            deallocate(omega_l)
            deallocate(l_core_x)
            deallocate(l_core_y)
            deallocate(l_core_z)
            
        end do
    
        write(*,*) 'Program complete.'
      
    end program main