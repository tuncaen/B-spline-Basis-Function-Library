module misc
    
    use Parameters, only:wp, eps
    use ogpf
    use vecfor

    implicit none
    
    interface find
        module procedure :: findindices
    end interface

    interface operator(.iseq.)
        module procedure real_equal
        module procedure char_equal
    end interface

    interface between
        module procedure real_between
        module procedure int_between
    end interface

    contains
    
    !# Plots 2D line graphs
    subroutine mlplot(x,xm,y,ym,ymg,fname,work,title,terminal,legend,showPlot,plotOpt,gpfname)
        real(wp), optional, intent(in) :: x(:)                      !< Horizontal axis data as vector
        real(wp), optional, intent(in) :: xm(:,:)                   !< Horizontal axis data as matrix
        real(wp), optional, intent(in) :: y(:)                      !< Vertical axis data as vector
        real(wp), optional, intent(in) :: ym(:,:)                   !< Vertical axis data as matrix
        real(wp), optional, intent(in) :: ymg(:,:,:)                 !< Vertical axis data as matrix
        character(*),   intent(in)   ,  optional    :: fname        !< Filename of the \*.png output, if *GNUPlot* terminal is `png`
        character(*),   intent(in)   ,  optional    :: work         !< Current work name to create folder
        character(*),   intent(in)   ,  optional    :: title        !< Title of the graph
        character(*),   intent(in)   ,  optional    :: terminal     !< *GNUPlot* terminal, default = ` wxt `
        character(*),   intent(in)   ,  optional    :: legend(:)
        character(*),   intent(in)   ,  optional    :: gpfname     !< *GNUPlot* script filename, input a name if you want to store script file
        logical, intent(in), optional       :: showPlot !< Key to open ` png ` files
        character(*),   intent(in)   ,  optional    :: plotOpt(:)   !< *GNUPlot* options
        integer :: i, j, k, iunit
        character(32), allocatable :: leg
        real(wp), allocatable :: xl(:,:)
        logical :: openfile
        logical :: legendon 
        character(200)  :: fn
        character(1024) :: ps
        
        type(gpf) :: gp
                    
        openfile = .false.; if(present(showPlot)) openfile = showPlot

        !set graph title
        if(present(title)) call gp%title(trim(title))
                
        !set gnuplot terminal
        ! if(present(terminal)) gp%term_type=trim(terminal) !etnc:gp%setterminal(terminal)
        
        !set output file name and fullpath
        fn = "Results"; if(present(fname)) fn = fname
        if(present(work)) then
            call gp%setoutput(folder="data/"//trim(work)//"/img/", fname=fn, term=terminal, gpfname=gpfname)
        else
            call gp%setoutput(folder="data/img/", fname=fn, term=terminal, gpfname=gpfname)
        end if

        ! Annotation: set title, xlabel, ylabel, line specification
        call gp%add_script('set xlabel "x"')
        call gp%add_script('set ylabel "y(x)"')
        !call gp%options('set format y "%.2sx10^{%S}"')
        !set plot options
        if(present(plotOpt))then
            do i=1,size(plotOpt)
                call gp%add_script(trim(plotOpt(i)))
            end do
        end if
        
        if(present(y))then
            if(present(x))then
                if(size(x,1)/=size(y,1)) then
                    print*, "error lplot: lenght of 'x' and 'y' variables does not match"
                    return
                end if
                call gp%plot(x,y,'with linespoints lt 2 pt 4')
            else
                print*, "error lplot: if dataset 'y' is a 1D array 'x' must be presented"
                return
            end if
        elseif(present(ym))then
            allocate(xl(size(ym,1),size(ym,2)))
            if(present(x))then
                if(size(x,1)/=size(ym,1)) then
                    print*, "error lplot: lenght of 'x' and 'ym' variables does not match"
                    return
                end if
                forall(i=1:size(xl,2)) xl(:,i) = x
            elseif(present(xm))then
                if(size(xm,1)/=size(ym,1) .or. size(xm,2)/=size(ym,2)) then
                    print*, "error lplot: lenght of 'xm' and 'ym' variables does not match"
                    return
                end if
                xl = xm
            end if
        elseif(present(ymg))then
            if(present(x))then
                if(size(x,1)/=size(ymg,1)) then
                    print*, "error lplot: lenght of 'x' and 'ymg' variables does not match"
                    return
                end if
            elseif(present(xm))then
                print*, "error lplot: to plot grouped dataset x-axis data should be the same"
                return
            end if
        endif
        
        if(present(ym) .and. (present(x) .or. present(xm)))then  ! if input has multiple data set

            if(present(gpfname))then
                open(newunit=iunit,file=trim(gp%folderpath)//trim(gpfname)//'_data.txt')
            else
                open(newunit=iunit,file='data.txt')
            end if
            legendon = .false.
            if(size(legend,1) ==  size(xl,2)) legendon = .true.
            do j=1,size(ym,2) ! do-for each set
                if(legendon)then
                    leg = trim(legend(j))
                    write(iunit,'(a,a,a)') '"',trim(leg),'"' 
                end if
                do i=1,size(xl,1)    !do for each data point
                    write(iunit,'(2(2x,e23.15e3))') xl(i,j), ym(i,j)
                end do
                write(iunit,*) new_line('a')
            end do
            
            !this script for exact solution
            !call gp%add_script('set linetype 1 lc rgb "black" lw 3 pt -1') 
            call gp%add_script('set key opaque') 
            !call gp%options('set offset graph 0.00, 0.00, 0.50, 0.50')
            
            if(present(gpfname))then
                write(ps,'("plot for [i=0:*] ",a,a,a," index i u 1:2 title columnheader(1)")') "'",trim(gpfname),"_data.txt'"
            else
                write(ps,'("plot for [i=0:*] ",a," index i u 1:2 title columnheader(1)")') "'data.txt'"
            end if
            call gp%add_script(ps)
            call gp%run_script()
            call gp%reset()

            close(iunit)

        elseif(present(ymg) .and. present(x))then  ! if input has grouped multiple data set

            if(present(gpfname))then
                open(newunit=iunit,file=trim(gp%folderpath)//trim(gpfname)//'_data.txt')
            else
                open(newunit=iunit,file='data.txt')
            end if
            legendon = .false.

            if(size(legend,1) ==  size(ymg,3)) legendon = .true. !legends are the groupnames

            do k=1,size(ymg,3)  ! do-for each group

                if(legendon)then
                    leg = trim(legend(k))
                    write(iunit,'(a,a,a)') '"',trim(leg),'"'
                end if
                ! do j=1,size(ymg,2) ! do-for each set

                    do i=1,size(x,1)    !do for each data point
                        write(iunit,'(20(2x,e23.15e3))') x(i), (ymg(i,j,k),j=1,size(ymg,2))
                    end do
                    write(iunit,*) new_line('a')

                ! end do

            end do
            
            !this script for exact solution
            !call gp%add_script('set linetype 1 lc rgb "black" lw 3 pt -1') 
            call gp%add_script('set key opaque')
            call gp%add_script('set colorsequence podo') 
            !call gp%options('set offset graph 0.00, 0.00, 0.50, 0.50')

            if(present(gpfname))then
                write(ps,'("plot for [i=0:*] ",a,a,a," index i u 1:2 title columnheader(1) lt i+1 lw 2 pt -1")') &
                "'",trim(gpfname),"_data.txt'"
            else
                write(ps,'("plot for [i=0:*] ",a," index i u 1:2 title columnheader(1) lt i+1 lw 2 pt -1")') "'data.txt'"
            end if
            ! ps = "plot for [i=0:*] 'data.txt' index i u 1:2 title columnheader(1) lt i+1 lw 2 pt -1"
            do k=2,size(ymg,2)
                write(ps,'(a,", for [i=0:*] ",a," index i u 1:",i0," notitle lt i+1 lw 2 pt -1")') trim(ps),"''",k+1
            end do

            ! "plot for [i=0:*] 'data.txt' index i u 1:2 title columnheader(1) ls 1 pt -1, '' u 1:3 notitle ls 1 pt -1"
            call gp%add_script(ps)
            call gp%run_script()
            call gp%reset()

            close(iunit)

        endif
        
        if(present(terminal) .and. openfile) then
            if(trim(terminal)=="png") then
                call execute_command_line(trim(gp%fullpath))
            end if
        end if
        
    end subroutine mlplot

    subroutine srfplot(x,y,z,xx,yy,zz,fname,work,title,terminal,showPlot,plotOpt,gpfname,type)
        real(wp), intent(in),optional :: x(:,:),y(:,:),z(:,:)
        real(wp), intent(in),optional :: xx(:,:,:),yy(:,:,:),zz(:,:,:)

        character(*),   intent(in)   ,  optional    :: fname        !< Filename of the \*.png output, if *GNUPlot* terminal is `png`
        character(*),   intent(in)   ,  optional    :: work         !< Current work name to create folder
        character(*),   intent(in)   ,  optional    :: title        !< Title of the graph
        character(*),   intent(in)   ,  optional    :: terminal     !< *GNUPlot* terminal, default = ` wxt `
        logical,        intent(in)   ,  optional    :: showPlot     !< Key to open ` png ` files
        character(*),   intent(in)   ,  optional    :: plotOpt(:)   !< *GNUPlot* options
        character(*),   intent(in)   ,  optional    :: gpfname     !< *GNUPlot* script filename, input a name if you want to store script file
        character(*),   intent(in)   ,  optional    :: type     !< *GNUPlot* script filename, input a name if you want to store script file
        !real(wp), allocatable :: xx(:,:),yy(:,:),zz(:,:)
        integer :: i, j, k, iunit
        character(200)  :: fn
        logical :: openfile, single, multi
        character(1024) :: ps
        
        type(gpf) :: gp
                 
        !single patch - multipatch plots
        single = present(x) .and. present(y) .and. present(z)
        multi = present(xx) .and. present(yy) .and. present(zz)

        !check sizes
        if(single)then
            if(.not.((size(x,1)==size(y,1)).and.(size(x,1)==size(z,1)))) return
            if(.not.((size(x,2)==size(y,2)).and.(size(x,2)==size(z,2)))) return
        elseif(multi)then
            if(.not.((size(xx,1)==size(yy,1)).and.(size(xx,1)==size(zz,1)))) return
            if(.not.((size(xx,2)==size(yy,2)).and.(size(xx,2)==size(zz,2)))) return
            if(.not.((size(xx,3)==size(yy,3)).and.(size(xx,3)==size(zz,3)))) return
        end if

        openfile = .false.; if(present(showPlot)) openfile = showPlot

        !set graph title
        if(present(title)) call gp%title(trim(title))
                
        !set gnuplot terminal
        ! if(present(terminal)) gp%term_type=trim(terminal)!etnc:gp%setterminal(terminal)
        
        !set output file name and fullpath 
        fn = "Results"; if(present(fname)) fn = fname
        if(present(work)) then
            call gp%setoutput(folder="data/"//trim(work)//"/img/", fname=fn, term=terminal, gpfname=gpfname)
        else
            call gp%setoutput(folder="data/img/", fname=fn, term=terminal, gpfname=gpfname)
        end if

        !set plot options
        if(present(plotOpt))then
            do i=1,size(plotOpt)
                call gp%add_script(trim(plotOpt(i)))
            end do
        end if
               
        if(present(gpfname))then
            open(newunit=iunit,file=trim(gp%folderpath)//trim(gpfname)//'_data.txt')
        else
            open(newunit=iunit,file='data.txt')
        end if

        if(single)then

            do j=1,size(x,2)
                do i=1,size(x,1)
                    write(iunit,'(3(e23.15e3,2x))') x(i,j), y(i,j), z(i,j)
                end do
                write(iunit,'(a)') 
            end do

        elseif(multi)then

            do k=1,size(xx,1)
                write(iunit,'(a,i0,a)') '"',k,'"'
                do j=1,size(xx,3)
                    do i=1,size(xx,2)
                        write(iunit,'(3(e23.15e3,2x))') xx(k,i,j), yy(k,i,j), zz(k,i,j)
                    end do
                    write(iunit,'(a)') 
                end do
                write(iunit,'(a)') 
            end do

        end if
        close(iunit)

        !graph title and output path
        !call gp%options('set offset graph 0.20, 0.20, 0.20, 0.20')
        !call gp%options('set dgrid3d 30,30,spline')
        !! call gp%options('set encoding utf8')
        ! call gp%add_script('set contour base')
        ! call gp%add_script('set cntrparam levels 14')
        if(type.iseq."contour")then
            call gp%options('set lmargin at screen 0.1')
            call gp%options('set rmargin at screen 0.8')     
            call gp%add_script('unset surface')     
            call gp%add_script('set view map')
        end if
        call gp%add_script(color_palettes('jet'))
        call gp%add_script('set pm3d') 
        ! call gp%add_script('set hidden3d')
        !call gp%add_script('splot "data.txt" using 1:2:3')
        
        if(present(gpfname))then
            write(ps,'("splot for [i=0:*] ",a,a,a," index i u 1:2:3 notitle")') "'",trim(gpfname),"_data.txt'"
        else
            write(ps,'("splot for [i=0:*] ",a," index i u 1:2:3 notitle")') "'data.txt'"
        end if

        call gp%add_script(ps)  
        !
        call gp%run_script()
        call gp%reset()

        !  call gp%contour(x, y, z, palette='jet' )
        
        if(present(terminal) .and. openfile) then
            if(trim(terminal)=="png") then
                call execute_command_line(trim(gp%fullpath))
            end if
        end if
         
    end subroutine srfplot

    ! subroutine surf4D(x,y,z,xx,yy,zz,fname,work,title,terminal,showPlot,plotOpt,gpfname)
    !     real(wp), intent(in),optional :: x(:,:),y(:,:),z(:,:)
    !     real(wp), intent(in),optional :: xx(:,:,:),yy(:,:,:),zz(:,:,:)

    !     character(*),   intent(in)   ,  optional    :: fname        !< Filename of the \*.png output, if *GNUPlot* terminal is `png`
    !     character(*),   intent(in)   ,  optional    :: work         !< Current work name to create folder
    !     character(*),   intent(in)   ,  optional    :: title        !< Title of the graph
    !     character(*),   intent(in)   ,  optional    :: terminal     !< *GNUPlot* terminal, default = ` wxt `
    !     logical,        intent(in)   ,  optional    :: showPlot     !< Key to open ` png ` files
    !     character(*),   intent(in)   ,  optional    :: plotOpt(:)   !< *GNUPlot* options
    !     character(*),   intent(in)   ,  optional    :: gpfname     !< *GNUPlot* script filename, input a name if you want to store script file
    !     !real(wp), allocatable :: xx(:,:),yy(:,:),zz(:,:)
    !     integer :: i, j, k, iunit
    !     character(200)  :: fn
    !     logical :: openfile, single, multi
    !     character(1024) :: ps
        
    !     type(gpf) :: gp
                 
    !     single = present(x) .and. present(y) .and. present(z)
    !     multi = present(xx) .and. present(yy) .and. present(zz)

    !     !check sizes
    !     if(single)then
    !         if(.not.((size(x,1)==size(y,1)).and.(size(x,1)==size(z,1)))) return
    !         if(.not.((size(x,2)==size(y,2)).and.(size(x,2)==size(z,2)))) return
    !     elseif(multi)then
    !         if(.not.((size(xx,1)==size(yy,1)).and.(size(xx,1)==size(zz,1)))) return
    !         if(.not.((size(xx,2)==size(yy,2)).and.(size(xx,2)==size(zz,2)))) return
    !         if(.not.((size(xx,3)==size(yy,3)).and.(size(xx,3)==size(zz,3)))) return
    !     end if

    !     openfile = .false.; if(present(showPlot)) openfile = showPlot

    !     !set graph title
    !     if(present(title)) call gp%title(trim(title))
                
    !     !set gnuplot terminal
    !     ! if(present(terminal)) gp%term_type=trim(terminal)!etnc:gp%setterminal(terminal)
        
    !     !set output file name and fullpath 
    !     fn = "Results"; if(present(fname)) fn = fname
    !     if(present(work)) then
    !         call gp%setoutput(folder="data/"//trim(work)//"/img/", fname=fn, term=terminal, gpfname=gpfname)
    !     else
    !         call gp%setoutput(folder="data/img/", fname=fn, term=terminal, gpfname=gpfname)
    !     end if

    !     !set plot options
    !     if(present(plotOpt))then
    !         do i=1,size(plotOpt)
    !             call gp%add_script(trim(plotOpt(i)))
    !         end do
    !     end if
               
    !     if(present(gpfname))then
    !         open(newunit=iunit,file=trim(gp%folderpath)//trim(gpfname)//'_data.txt')
    !     else
    !         open(newunit=iunit,file='data.txt')
    !     end if

    !     if(single)then

    !         do j=1,size(x,2)
    !             do i=1,size(x,1)
    !                 write(iunit,'(3(e23.15e3,2x))') x(i,j), y(i,j), z(i,j)
    !             end do
    !             write(iunit,'(a)') 
    !         end do

    !     elseif(multi)then

    !         do k=1,size(xx,1)
    !             write(iunit,'(a,i0,a)') '"',k,'"'
    !             do j=1,size(xx,3)
    !                 do i=1,size(xx,2)
    !                     write(iunit,'(3(e23.15e3,2x))') xx(k,i,j), yy(k,i,j), zz(k,i,j)
    !                 end do
    !                 write(iunit,'(a)') 
    !             end do
    !             write(iunit,'(a)') 
    !         end do

    !     end if
    !     close(iunit)

    !     !graph title and output path
    !     !call gp%options('set offset graph 0.20, 0.20, 0.20, 0.20')
    !     !call gp%options('set dgrid3d 30,30,spline')
    !     !! call gp%options('set encoding utf8')
    !     call gp%options('set lmargin at screen 0.1')
    !     call gp%options('set rmargin at screen 0.8')     
    !     ! call gp%add_script('set contour base')
    !     ! call gp%add_script('set cntrparam levels 14')
    !     call gp%add_script('unset surface')     
    !     call gp%add_script('set view map')
    !     !call gp%add_script('set hidden3d')
    !     call gp%add_script(color_palettes('jet'))
    !     call gp%add_script('set pm3d') 
    !     !call gp%add_script('splot "data.txt" using 1:2:3')
        
    !     if(present(gpfname))then
    !         write(ps,'("splot for [i=0:*] ",a,a,a," index i u 1:2:3 notitle")') "'",trim(gpfname),"_data.txt'"
    !     else
    !         write(ps,'("splot for [i=0:*] ",a," index i u 1:2:3 notitle")') "'data.txt'"
    !     end if

    !     call gp%add_script(ps)  
    !     !
    !     call gp%run_script()
    !     call gp%reset()

    !     !  call gp%contour(x, y, z, palette='jet' )
        
    !     if(present(terminal) .and. openfile) then
    !         if(trim(terminal)=="png") then
    !             call execute_command_line(trim(gp%fullpath))
    !         end if
    !     end if
        
    ! end subroutine surf4D


    ! This routine finds the indices of .true. value in the logical input
    ! etnc: add usage example    
    pure function findindices(n,tf) result(indices)

        ! inlet variables
        integer,intent(in):: n      ! dimension of logical vector
        logical,intent(in):: tf(n)  ! logical vector (true or false)
        ! outlet variables
        integer npos                ! number of "true" conditions
        integer pos(n)              ! position of "true" conditions
        ! internal variables
        integer i                   ! counter
        integer v(n)                ! vector of all positions
        integer, dimension(:), allocatable :: indices

        pos = 0                     ! initialize pos
        forall(i=1:n)   v(i) = i    ! enumerate all positions
        npos  = count(tf)           ! count the elements of tf that are .true.
        pos(1:npos)= pack(v, tf)    ! with pack function, verify position of true conditions

        allocate(indices(npos))
        indices = pos(1:npos)

    end function findindices

    pure function lowcase(s) result(t)
        ! Returns string 's' in lowercase
        character(*), intent(in) :: s
        character(len(s)) :: t
        integer :: i, diff
        t = s; diff = ichar('A')-ichar('a')
        do i = 1, len(t)
            if (ichar(t(i:i)) >= ichar('A') .and. ichar(t(i:i)) <= ichar('Z')) then
                ! if uppercase, make lowercase
                t(i:i) = char(ichar(t(i:i)) - diff)
            end if
        end do
    end function lowcase

    !# logical function to compare real values
    pure elemental logical function real_equal(r1,r2)
        real(wp), intent(in) :: r1
        real(wp), intent(in) :: r2
        real(wp), parameter :: eps_wp  = epsilon(eps_wp)
        real(wp), parameter :: eps_wp3 = 3.0_wp * epsilon(eps_wp)

        !real_equal = abs(r1-r2) < eps
        real_equal = abs(r1-r2) <= max( abs(r1), abs(r2) ) * eps_wp3

    end function real_equal
    
    !# logical function to compare character keys
    pure elemental logical function char_equal(c1,c2)
        character(*), intent(in) :: c1
        character(*), intent(in) :: c2

        char_equal = lowcase(trim(c1))==lowcase(trim(c2))

    end function char_equal
    
    pure elemental logical function isClose(r1,r2,tol)
        real(wp), intent(in) :: r1
        real(wp), intent(in) :: r2
        real(wp), intent(in), optional :: tol
        real(wp) :: toler
        ! real(wp), parameter :: eps_wp  = epsilon(eps_wp)

        toler = 1e-8_wp; if(present(tol)) toler=tol
        isClose = abs(r1-r2) <=  toler !max( abs(r1), abs(r2) ) * 3.0_wp *

    end function isClose

    !# Prints runtime
    subroutine run_time(time0,time1,str)

        implicit none
        real(wp),intent(in) :: time0,time1
        
        ! local
        integer :: sec,min,hr
        character(len=*),intent(in) :: str
        sec = int(time1-time0) ; hr = sec/3600 ; sec = sec-hr*3600 ; min = sec/60 ; sec = sec-min*60
        write (*,'(a,2(i2.2,":"),i2.2,/)') "runtime for "//trim(str)//": ",hr,min,sec   
    
    end subroutine run_time

    !# <a href="https://stackoverflow.com/questions/58938347/how-do-i-replace-a-character-in-the-string-with-another-charater-in-fortran" 
    !   target="_blank">Find and replace characters in Fortran strings</a> 
    pure recursive function replaceStr(string,search,substitute) result(modifiedString)
        character(len=*), intent(in)  :: string, search, substitute
        character(len=:), allocatable :: modifiedString
        integer                       :: i, stringLen, searchLen
        stringLen = len(string)
        searchLen = len(search)
        if (stringLen==0 .or. searchLen==0) then
            modifiedString = ""
            return
        elseif (stringLen<searchLen) then
            modifiedString = string
            return
        end if
        i = 1
        do
            if (string(i:i+searchLen-1)==search) then
                modifiedString = string(1:i-1) // substitute // replaceStr(string(i+searchLen:stringLen),search,substitute)
                exit
            end if
            if (i+searchLen>stringLen) then
                modifiedString = string
                exit
            end if
            i = i + 1
            cycle
        end do
    end function replaceStr

    pure function complement(a,b) result(r)
        integer, intent(in) :: a(:)
        integer, intent(in) :: b(:)
        integer, allocatable :: r(:)
        integer :: i,j,n,m

        n=size(a)
        m=size(b)
        allocate(r(n-m)); r=0

        j = 0
        do i=1,n
            if(any(b == a(i)))cycle
            j=j+1
            if(j>n-m) error stop "complement error"
            r(j) = a(i) 
        end do        
        
    end function complement

    pure function union(a, b) result(u)
        integer, intent(in) :: a(:), b(:)
        integer, allocatable :: u(:)
        integer :: i, k, m, n
        integer :: indices(size(b))

        indices = 0
        m = size(a); n = size(b)

        k = 0
        do i = 1, n
            if(any(a==b(i))) cycle
            k = k + 1
            indices(k) = i
        end do

        allocate(u(m+k))
        u(1:m) = a
        u(m+1:m+k) = b(indices(1:k))
        
    end function union

    pure function real_between(x,a,b,ex) result(r)
        real(wp), intent(in) :: x, a, b
        logical, intent(in), optional :: ex
        logical :: r
        logical :: exclude

        exclude = .false.; if(present(ex)) exclude = ex

        if(exclude)then
            r = (x>a) .and. (x<b)
        else
            r = (x>=a) .and. (x<=b)
        end if

    end function real_between

    pure function int_between(x,a,b,ex) result(r)
        integer, intent(in) :: x, a, b
        logical, intent(in), optional :: ex
        logical :: r
        logical :: exclude

        exclude = .false.; if(present(ex)) exclude = ex

        if(exclude)then
            r = (x>a) .and. (x<b)
        else
            r = (x>=a) .and. (x<=b)
        end if

    end function int_between

    pure function toVector(v) result(r)
        real(wp), intent(in) :: v(:)
        type(vector) :: r
        integer :: n

        n=size(v)
        if(n==1)then
            r = v(1)*ex
        elseif(n==2)then
            r = v(1)*ex + v(2)*ey
        elseif(n>=3)then
            r = v(1)*ex + v(2)*ey + v(3)*ez
        end if

    end function toVector

    ! subroutine meshgrid(x,y,xgv,ygv, ierr)
    !     !..............................................................................
    !     !meshgrid generate mesh grid over a rectangular domain of [xmin xmax, ymin, ymax]
    !     ! Inputs:
    !     !     xgv, ygv are grid vectors in form of full grid data
    !     ! Outputs:
    !     !     X and Y are matrix each of size [ny by nx] contains the grid data.
    !     !     The coordinates of point (i,j) is [X(i,j), Y(i,j)]
    !     !     ierr: The error flag
    !     !     """
    !     !     # Example
    !     !     # call meshgrid(X, Y, [0.,1.,2.,3.],[5.,6.,7.,8.])
    !     !     # X
    !     !     # [0.0, 1.0, 2.0, 3.0,
    !     !     #  0.0, 1.0, 2.0, 3.0,
    !     !     #  0.0, 1.0, 2.0, 3.0,
    !     !     #  0.0, 1.0, 2.0, 3.0]
    !     !     #
    !     !     #Y
    !     !     #[ 5.0, 5.0, 5.0, 5.0,
    !     !     #  6.0, 6.0, 6.0, 6.0,
    !     !     #  7.0, 7.0, 7.0, 7.0,
    !     !     #  8.0, 8.0, 8.0, 8.0]
    !     !..............................................................................
    !     ! Rev 0.2, Feb 2018
    !     ! New feature added: xgv and ygv as full grid vector are accepted now

    !     ! Arguments
    !     real(wp), intent(out), allocatable  :: x(:,:)
    !     real(wp), intent(out), allocatable  :: y(:,:)
    !     real(wp), intent(in)                :: xgv(:) ! x grid vector [start, stop, step] or [start, stop]
    !     real(wp), intent(in),  optional     :: ygv(:) ! y grid vector [start, stop, step] or [start, stop]
    !     integer,  intent(out), optional     :: ierr   ! the error value

    !     ! Local variables
    !     integer:: sv
    !     integer:: nx
    !     integer:: ny
    !     logical:: only_xgv_available

    !     ! Initial setting
    !     only_xgv_available  = .false.
    !     sv=0 !Assume no error

    !     nx=size(xgv, dim=1)

    !     if (present(ygv)) then
    !         ny = size(ygv, dim=1)
    !     else
    !         only_xgv_available=.true.
    !         ny=nx
    !     end if

    !     allocate(x(ny,nx),y(ny,nx),stat=sv)
    !     if (sv /=0) then
    !         print*, "allocataion erro in meshgrid"
    !         stop
    !     end if

    !     x(1,:)    = xgv
    !     x(2:ny,:) = spread(xgv, dim=1, ncopies=ny-1)

    !     if (only_xgv_available) then
    !         y=transpose(x)
    !     else
    !         y(:,1)    = ygv
    !         y(:,2:nx) = spread(ygv,dim=2,ncopies=nx-1)
    !     end if

    !     if (present(ierr)) then
    !         ierr=sv
    !     end if

    ! end subroutine meshgrid

end module misc