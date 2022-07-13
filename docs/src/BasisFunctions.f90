module BasisFunctions

    !< The module defining spline basis objects and their arguments

    use Parameters
    use ogpf
    use misc
    
    implicit none
        
    private

    public :: basis

    !> **Spline basis function class**  
    type :: basis
        !<  
        integer                 ::  p       !< Polynomial degree \(p\)
        integer                 ::  nk      !< Number of the knots
        integer                 ::  m       !< The highest index of the knots ` ( nk-1 ) `
        integer                 ::  nb      !< Number of basis functions
        integer                 ::  n       !< The highest index of the basis functions ` ( nb-1 ) `
        real(wp),   allocatable ::  kv(:)   !< [Open knot vector](../page/01.fooiga/eqns.html#openvec)  
                                            !< ` kv(0:m) ` \( = \left\{u_{0},..., u_{m}\right\} \)
        real(wp),   allocatable ::  unik(:) !< Unique knot values
        integer,    allocatable ::  mult(:) !< Multiplicity of each unique knot value
        integer                 ::  minreg  !< Minimum reqularity of the basis functions,` = me%p-max(me%mult(2:ubound(me%mult)-1)  `
        integer                 ::  maxreg  !< Maximum reqularity of the basis functions,` = me%p-min(me%mult(2:ubound(me%mult)-1)  `
        integer                 ::  nbez    !< Number of elements with nonzero measure, i.e., number of Bezier segments
    contains
        !< Procedures assume basis with **open knot vector**  
        procedure, public   :: FindSpan                 !< **Knot span index search                                                                     **
        procedure, public   :: FindMult                 !< **Knot multiplicity computation                                                              **
        procedure, public   :: FindSpanMult             !< **Knot span and multiplicity computation                                                     **
        procedure, public   :: BasisFuns                !< **Computation of the non-vanishing basis functions                                           **
        procedure, public   :: OneBasisFun              !< **Computation of a single basis function                                                     **
        procedure, public   :: AllBasisFuns             !< **Computation of a all non-zero basis function of all degrees from \(0\) up to \(p\)         **
        procedure, public   :: DersBasisFuns            !< **Computation of the derivatives of the non-vanishing basis functions                        **
        procedure, public   :: DersOneBasisFun          !< **Computation of the the derivatives of a single basis function                              **
        procedure, public   :: InsertKnot               !< **Knot insertion of a given knot into the knot vector                                        **
        procedure, public   :: UpdateBasis              !< **Updating the basis (the same contents with constructor)                                    **
        procedure, public   :: BSegments                !< **Computation of the segments' data of the basis                                             **
        procedure, public   :: ParLen                   !< **Get the length of parameter space                                                          **
        procedure, public   :: GrevilleAbscissae        !< **Computation of the Greville abscissae                                                      **  
        procedure, public   :: Plot => PlotBasisFunc    !< **Plot basis functions                                                                       **
        procedure, public   :: get_nquad                !< **Get the number of quadrature points                                                        **
    end type basis

    
    !***************************************************************************************************************************
    ! **OBJECT CONSTRUCTOR**  
    ! - Basis       (Int)  
    ! =>
        !> **Spline basis constructor**  
        interface basis
            procedure basisConstructor
        end interface

    !***************************************************************************************************************************

    contains


    !***************************************************************************************************************************
    ! **METHODS OF `BASIS` TYPE**  
    ! - PlotBasisFunc       (Sub)  
    ! - FindSpan            (Fn)  
    ! - FindMult            (Fn)  
    ! - FindSpanMult        (Sub)  
    ! - BasisFuns           (Sub)  
    ! - OneBasisFun         (Fn)    
    ! - AllBasisFuns        (Sub)  
    ! - DersBasisFuns       (Sub)  
    ! - DersOneBasisFun     (Sub)  
    ! - InsertKnot          (Sub)  
    ! - UpdateBasis         (Sub)  
    ! - BSegments           (Sub)  
    ! - ParLen              (Fn)  
    ! - GrevilleAbscissae   (Fn)  
    ! - get_nquad           (Fn)  
    ! =>

        !# {!module(basisfunctions)/type(basis)/PlotBasisFunc.md!}
        subroutine PlotBasisFunc(me,first,last,d,plotRes,title,work,fname,terminal,showPlot)
            ! Variables
            class(basis), intent(in)            :: me
            integer, intent(in), optional       :: first    !< First index to be plotted
            integer, intent(in), optional       :: last     !< Last index to be plotted
            integer, intent(in), optional       :: d        !< Derivatives degree
            integer, intent(in), optional       :: plotRes  !< Resolution of the calculation, Default = 100
            character(*), intent(in), optional  :: title    !< Title of the graph, Default = "Basis functions"
            character(*), intent(in), optional  :: work     !< Current work name to create folder
            character(*), intent(in), optional  :: fname    !< Filename of the \*.png output, if *GNUPlot* terminal is `png`. Default = "BasisFunctions"
            character(*), intent(in), optional  :: terminal !< *GNUPlot* terminal, Default = ` wxt `
            logical, intent(in), optional       :: showPlot !< Key to open created ` png ` files
            !locals
            integer :: i, j, k, m, r, f, l, dd, iunit
            character(200) :: fn, ptitle, fmt1, fmt2, yaxis
            character(1024) :: pl 
            real(wp) ::du, u
            real(wp), allocatable :: array(:), bf(:,:)
            logical :: openfile
            type(gpf):: gp
            
            !set gnuplot terminal
            !if(present(terminal)) gp%term_type=trim(terminal)
            openfile    = .false.   ; if(present(showPlot))     openfile    = showPlot
            r = 100; if(present(plotRes)) r=plotRes
            f = 0; if(present(first)) f=first
            l = me%n; if(present(last)) l=last
            dd = 0; if(present(d)) dd=d
            allocate(bf(0:me%p,0:dd))
            allocate(array(0:l-f))

            write(fn, '("_p",i0,"_d",i0,"_idx",i0,"_to_",i0)') me%p,dd,f,l
            if(present(fname)) then
                fn = fname // trim(fn)
            else
                fn = "BasisFunctions" // trim(fn)
            end if
 
            if(present(work)) then
                call gp%setoutput(folder="data/"//trim(work)//"/img/", fname=fn, term=terminal)
            else
                call gp%setoutput(folder="data/img/", fname=fn, term=terminal)
            end if
            
            !set plot title
            ptitle = "Basis functions"; if(present(title)) ptitle=title
            call gp%title(trim(ptitle))
            
            !set headers format
            write(fmt1,'(a,i2.1,a)') '("x-values",', l-f+1, '(1x,"N_{",i0,"}"))'
            !set data format
            write(fmt2,'(a,i2.1,a)') '(', l-f+2, '(2x,e23.15E3))'

            open(newunit=iunit,file='data.txt')
            !write headers
            write(iunit,trim(fmt1)) (j,j=f,l) 
            
            !calculation and writing data
            du = ( me%kv(me%m) - me%kv(0) ) / r
            do i=0,r
                u = i * du
                k = me%FindSpan(u)
                call me%DersBasisFuns(u,dd,bf,k)
                array = 0.0_wp; m = 0
                do j=f,l
                    if(k-me%p <= j .and. j <= k)then
                        array(j) = bf(m,dd)
                        m = m+1
                    end if
                end do
                write(iunit,fmt2) u, (array(j),j=f,l)
            end do
            
            call gp%xlabel('u')
            if(dd>0)then
                write(yaxis,'("N_{i}^{(",i0,")}(u)")') dd
            else
                yaxis = "N_{i}(u)"
            end if
            call gp%add_script('set ylabel "'//trim(yaxis)//'" rotate by 0')


            call gp%add_script('set tics') 
            call gp%add_script('set colorsequence podo') 
            
            !set plotting script
            write(pl,'("plot for [i=2:",i0,"] ",a," using 1:i w lines t columnheader(i) lw 2.5")') l-f+2, '"data.txt"'
            call gp%add_script(pl)
            call gp%run_script()
            call gp%reset()
            if(present(terminal) .and. openfile) then
                if(trim(terminal)=="png") then
                    call execute_command_line(trim(gp%fullpath))
                end if
            end if
            close(iunit,status='delete')
            
        end subroutine PlotBasisFunc
        ! .......................................................

        !# {!module(basisfunctions)/type(basis)/FindSpan.md!}
        pure elemental function FindSpan(me,u) result (span)
            class(basis), intent(in)    :: me
            real(wp), intent(in)        :: u        !< Given \( u \) value
            integer                     :: span     !< Knot span index of \( u \) value
            !locals
            integer :: low, high    !< indices for binary search
            
            associate( p => me%p, Ui => me%kv, n => me%n ) !get rid of % signs
            !/* Special cases */
            if (u >= Ui(n+1)) then
                span = me%n
                return
            end if
            if (u <= Ui(p)) then
            span = p
            return
            end if
            !/* Do binary search */
            low = p; high = n+1
            span = (low + high) / 2
            do while (u < Ui(span) .or. u >= Ui(span+1))
                if (u < Ui(span)) then
                    high = span
                else
                    low  = span
                end if
                span = (low + high) / 2
            end do
            end associate
            
        end function FindSpan
        ! .......................................................

        !# {!module(basisfunctions)/type(basis)/FindMult.md!}
        pure elemental function FindMult(me,i,u) result (mult)
            class(basis), intent(in)    :: me
            real(wp), intent(in)        :: u        !< Given \( u \) value
            integer, intent(in)         :: i        !< Starting index for search 
            integer                     :: mult     !< Multiplicity of \( u \) value
            !locals
            integer :: j
            
            associate( p => me%p, Ui => me%kv )
                mult = 0
                do j = -p, p+1
                    if (u == Ui(i+j)) mult = mult + 1
                end do
            end associate

        end function FindMult
        ! .......................................................
            
        !# {!module(basisfunctions)/type(basis)/FindSpanMult.md!}
        pure elemental subroutine FindSpanMult(me,u,k,s)
            class(basis), intent(in)    :: me
            real(wp), intent(in)        :: u        !< Given \( u \) value
            integer, intent(out)        :: k        !< Knot span index of \( u \) value
            integer, intent(out)        :: s        !< Multiplicity of \( u \) value
            
            k = me%FindSpan(u)
            s = me%FindMult(k,u)    
        end subroutine FindSpanMult
        ! .......................................................

        !# {!module(basisfunctions)/type(basis)/BasisFuns.md!}
        pure subroutine BasisFuns(me,u,N,span)
            class(basis), intent(in)                :: me
            real(wp),       intent(in)              :: u            !< Given \( u \) value
            integer,        intent(in), optional    :: span         !< Knot span where \( u \) lies on
            real(wp),       intent(out)             :: N(0:me%p)    !< \( =  \{ N_{i-p,p}(u) \ldots N_{i,p}(u)\} \)
            !locals
            integer     :: i    !< Knot span index
            integer     :: j, r
            real(wp)    :: left(me%p), right(me%p), saved, temp!, w
            
            associate( p => me%p, Ui => me%kv ) !get rid of % signs
                
            ! set the knot span
            if(present(span))then
                i = span
            else
                i = me%FindSpan(u)
            end if

            N = 0.0_wp; right=0.0_wp; left=0.0_wp
            
            N(0) = 1.0_wp
            do j = 1, p
            left(j)  = u - Ui(i+1-j)
            right(j) = Ui(i+j) - u
            saved = 0.0_wp
            do r = 0, j-1
                temp = N(r) / (right(r+1) + left(j-r))
                N(r) = saved + right(r+1) * temp
                saved = left(j-r) * temp
            end do
            N(j) = saved
            end do

            ! if(me%rational)then
            !     w = sum(N(:)*me%w(i-p:i))
            !     do j = 0, me%p
            !         N(j) = N(j)*me%w(j) / w
            !     end do
            ! end if
                
            end associate
            
        end subroutine BasisFuns
        ! .......................................................

        !#  {!module(basisfunctions)/type(basis)/OneBasisFun.md!}
        pure elemental function OneBasisFun(me,i,u) result(Nip)
            class(basis), intent(in)    :: me
            integer,        intent(in)  :: i    !< Basis function index   
            real(wp),       intent(in)  :: u    !< Given \( u \) value
            real(wp)                    :: Nip  !< \( N_{i,p}(u)\)
            !locals
            integer     :: j, k
            real(wp)    :: N(0:me%p), Uleft, Uright, saved, temp
            !initialize arrays

            ! if(me%rational)then
            !     print*, "'OneBasisFun' does not support for rational basis, you already have to compute"
            !     print*, "       all non-zero function for rational basis. Call 'BasisFuns'"
            !     Nip = 0.0_wp
            !     return
            ! end if
            
            N = 0.0_wp
            associate( p => me%p, Ui => me%kv, m => me%m ) !get rid of % signs
            
                if (( i == 0 .and. u == Ui(0)) .or. (i == m-p-1 .and. u == Ui(m))) then !/* Special cases */
                    Nip = 1.0_wp; return
                end if
                
                if (u < Ui(i) .or. u >= Ui(i+p+1)) then ! /* Local property */
                    Nip = 0.0_wp; return
                end if
                
                do j=0,p ! /* Initialize zeroth-degree functs */
                    if (u >= Ui(i+j) .and. u < Ui(i+j+1))then
                        N(j) = 1.0_wp
                    else 
                        N(j) = 0.0_wp
                    end if
                end do
                
                do k=1,p ! /* Compute triangular table */
                    if (N(0) == 0.0_wp)then
                        saved = 0.0_wp
                    else 
                        saved = (( u-Ui(i) ) * N(0) ) / ( Ui(i+k) - Ui(i) )
                    end if
                    do j=0,p-k
                        Uleft = Ui(i+j+1)
                        Uright = Ui(i+j+k+1)
                        if (N(j+1) == 0.0_wp)then
                            N(j) = saved 
                            saved = 0.0_wp
                        else
                            temp = N(j+1)/(Uright-Uleft)
                            N(j) = saved+(Uright-u)*temp
                            saved = (u-Uleft)*temp
                        end if
                    end do
                end do
                Nip = N(0)
            end associate
            
        end function OneBasisFun
        ! .......................................................

        !# {!module(basisfunctions)/type(basis)/AllBasisFuns.md!}
        pure subroutine AllBasisFuns(me,u,N,span)
            class(basis), intent(in)            :: me
            real(wp),   intent(in)              :: u        !< Given \( u \) value
            integer,    intent(in), optional    :: span     !< Knot span where \( u \) lies on
            real(wp),   intent(out)             :: N(0:me%p, 0:me%p)   !<
                                                !< $$ = \begin{bmatrix} {N_{i - p,0}(u)} & \cdots &{N_{i - p,p}(u)}\\
                                                !< \vdots & \ddots & \vdots \\ {N_{i,0}(u)} & \cdots &{N_{i,p}(u)} \end{bmatrix} $$
            !locals
            integer     :: i    !< Knot span index
            integer     :: j, k, r
            real(wp)    :: left(me%p), right(me%p), saved, temp
            

            associate( p => me%p, Ui => me%kv ) !get rid of % signs
                
            ! set the knot span
            if(present(span))then
                i = span
            else
                i = me%FindSpan(u)
            end if

            N = 0.0_wp; right=0.0_wp; left=0.0_wp
            
            N(0,0) = 1.0_wp
            do k = 1, p
            do j = 1, k
            left(j)  = u - Ui(i+1-j)
            right(j) = Ui(i+j) - u
            saved = 0.0_wp
            do r = 0, j-1
                temp = N(r,k-1) / (right(r+1) + left(j-r))
                N(r,k) = saved + right(r+1) * temp
                saved = left(j-r) * temp
            end do
            N(j,k) = saved
            end do
            end do
                
            end associate

        end subroutine AllBasisFuns
        ! .......................................................

        !# {!module(basisfunctions)/type(basis)/DersBasisFuns.md!}
        pure subroutine DersBasisFuns(me,u,n,ders,span)
            class(basis), intent(in)            :: me
            real(wp),   intent(in)              :: u            !< Given \( u \) value
            integer,    intent(in)              :: n            !< Number of derivatives (\(n \leq p\))
            integer,    intent(in), optional    :: span         !< Knot span where \( u \) lies on
            real(wp),   intent(out)             :: ders(0:me%p, 0:n)    !<
                                                !< $$ = \begin{bmatrix} {N_{i - p,p}^{(0)}(u)} & \cdots &{N_{i - p,p}^{(k)}(u)}\\
                                                !< \vdots & \ddots & \vdots \\ {N_{i,p}^{(0)}(u)} & \cdots &{N_{i,p}^{(k)}(u)} \end{bmatrix} $$  
            !locals
            integer :: i    !< Knot span index
            integer :: j, k, r, s1, s2, rk, pk, j1, j2
            real(wp) :: saved, temp, d
            real(wp) :: left(me%p), right(me%p)
            real(wp) :: ndu(0:me%p,0:me%p), a(0:1,0:me%p)
            
            ders = 0.0_wp
            associate( p => me%p, Ui => me%kv ) !get rid of % signs
            ! set the knot span
            if(present(span))then
                i = span
            else
                i = me%FindSpan(u)
            end if
                
            ndu(0,0) = 1.0_wp
            do j = 1, p
                left(j)  = u - Ui(i+1-j)
                right(j) = Ui(i+j) - u
                saved = 0.0_wp
                do r = 0, j-1
                    !/* Lower triangle */
                    ndu(j,r) = right(r+1) + left(j-r) 
                    temp = ndu(r,j-1) / ndu(j,r)
                    !/* Upper triangle */
                    ndu(r,j) = saved + right(r+1) * temp
                    saved = left(j-r) * temp
                end do
                ndu(j,j) = saved
            end do
            
            !/* Load the basis functions */
            ders(:,0) = ndu(:,p)
            
            !/* This section computes the derivatives (Eq. [2.9]) */
            do r = 0, p !/* Loop over function index */
                s1 = 0; s2 = 1 !/* Alternate rows in array a */
                a(0,0) = 1.0_wp
                !/* Loop to compute kth derivative */
                do k = 1, n
                    d = 0.0_wp
                    rk = r-k; pk = p-k;
                    if (r >= k) then
                        a(s2,0) = a(s1,0) / ndu(pk+1,rk)
                        d =  a(s2,0) * ndu(rk,pk)
                    end if
                    if (rk > -1) then
                        j1 = 1
                    else
                        j1 = -rk
                    end if
                    if (r-1 <= pk) then
                        j2 = k-1
                    else
                        j2 = p-r
                    end if
                    do j = j1, j2
                        a(s2,j) = (a(s1,j) - a(s1,j-1)) / ndu(pk+1,rk+j)
                        d =  d + a(s2,j) * ndu(rk+j,pk)
                    end do
                    if (r <= pk) then
                        a(s2,k) = - a(s1,k-1) / ndu(pk+1,r)
                        d =  d + a(s2,k) * ndu(r,pk)
                    end if
                    ders(r,k) = d
                    j = s1; s1 = s2; s2 = j !/* Switch rows */
                end do
            end do
            !/* Multiply through by the correct factors */
            !/* (Eq. [2.9]) */
            r = p
            do k = 1, n
                ders(:,k) = ders(:,k) * r
                r = r * (p-k)
            end do
            end associate

        end subroutine DersBasisFuns
        ! .......................................................

        !# {!module(basisfunctions)/type(basis)/DersOneBasisFun.md!}
        pure subroutine DersOneBasisFun(me,i,u,n,ders)
            class(basis),   intent(in)  :: me
            integer,        intent(in)  :: i            !< Basis function index   
            real(wp),       intent(in)  :: u            !< Given \( u \) value
            integer,        intent(in)  :: n            !< Number of derivatives (\(n \leq p\))
            real(wp),       intent(out) :: ders(0:n)    !< \( = \{ N_{i,p}^{(0)}(u), \ldots, N_{i,p}^{(k)}(u) \} \)
            !locals
            integer                 :: m
            integer                 :: j, jj, k
            real(wp)                :: ND(0:me%p), Ni(0:me%p,0:n), Uleft, Uright, saved, temp
            
            
            associate( p => me%p, Ui => me%kv ) !get rid of % signs 
            
            ! set highest knot vector index 
            m = me%nk-1
            ND = 0.0_wp
            Ni = 0.0_wp
            ders = 0.0_wp
            
            ! if knot is outside of span range
            if (u < Ui(i) .or. u >= Ui(i+p+1))then ! /* Local property */
                do k=0,n 
                    ders(k) = 0.0_wp
                end do
                return
            end if
            
            do j=0,p ! /* Initialize zeroth-degree functions */
                if (u >= Ui(i+j) .and. u < Ui(i+j+1)) then
                    Ni(j,0) = 1.0_wp
                else 
                    Ni(j,0) = 0.0_wp
                end if
            end do
            
            do k=1,p ! / * Compute full triangular table * /
                ! Detecting zeros saves computations
                if (Ni(0,k-1) == 0.0_wp)then
                    saved = 0.0_wp    
                else 
                    saved = ((u-Ui(i)) * Ni(0,k-1)) / (Ui(i+k)-Ui(i))
                end if
                do j=0,p-k
                    Uleft = Ui(i+j+1)
                    Uright = Ui(i+j+k+1)
                    ! Zero detection
                    if (Ni(j+1,k-1) == 0.0_wp)then
                        Ni(j,k) = saved 
                        saved = 0.0_wp
                    else
                        temp = Ni(j+1,k-1)/(Uright-Uleft)
                        Ni(j,k) = saved+(Uright-u)*temp
                        saved = (u-Uleft)*temp
                    end if
                end do
            end do
            
            ders(0) = Ni(0,p) ! /* The function value */
            
            do k=1,n ! /* Compute the derivatives */
                ND = 0.0_wp
                do j=0,k ! /* Load appropriate column */
                    ND(j) = Ni(j,p-k) 
                end do
                do jj=1,k ! /* Compute table of width k */
                    if (ND(0) == 0.0_wp)then
                        saved = 0.0_wp
                    else
                        saved = ND(0) / (Ui(i+p-k+jj) - Ui(i))
                    end if
                    do j=0,k-jj
                        Uleft = Ui(i+j+1)
                        Uright = Ui(i+j+p-k+jj+1) !Wrong in The NURBS Book: -k is missing.
                        if (ND(j+1) == 0.0_wp)then
                            ND(j) = (p-k+jj)*saved
                            saved = 0.0_wp
                        else
                            temp = ND(j+1)/(Uright-Uleft);
                            ND(j) = (p-k+jj)*(saved-temp);
                            saved = temp;
                        end if
                    end do
                end do
                ders(k) = ND(0) ! /* kth derivative */
            end do
                
            end associate
            
        end subroutine DersOneBasisFun
        ! .......................................................

        !# {!module(basisfunctions)/type(basis)/InsertKnot.md!}
        pure subroutine InsertKnot(me,u,r,span)
            class(basis), intent(inout)     :: me
            real(wp), intent(in)            :: u    !< knot value \(\bar{u}\)
            integer, intent(in)             :: r    !< multiplicity of new knot value
            integer, intent(in), optional   :: span !< knot span
            !locals
            real(wp) :: kvtemp(0:me%m)
            integer  :: k
            
            if(present(span))then
                k = span
            else
                k = me%FindSpan(u)
            end if
            
            kvtemp = me%kv
            deallocate(me%kv)
            allocate(me%kv(0:me%m+r))
            me%kv(0:k) = kvtemp(0:k)
            me%kv(k+1:k+r) = u
            me%kv(k+r+1:me%m+r) = kvtemp(k+1:me%m)
            
            !update basis properties
            me%nk = me%nk + r
            me%m = me%m + r
            me%nb = me%nb + r
            me%n = me%n + r
        end subroutine InsertKnot
        ! .......................................................

        !# {!module(basisfunctions)/type(basis)/UpdateBasis.md!}
        pure subroutine UpdateBasis(me,nk,k,p)
            class(basis), intent(inout)     :: me
            integer, intent(in)             :: nk       !< Number of the knots
            real(wp), intent(in)            :: k(nk)    !< New knot vector 
            integer, intent(in), optional   :: p        !< New basis degree
            !locals
            ! real(wp), allocatable       :: uniq(:)
            integer                     :: i

            deallocate(me%kv)
            if(present(p)) me%p = p
            me%nk   = nk
            me%m    = nk - 1
            me%nb   = nk - ( me%p + 1 )     !number of basis functions or control points
            me%n    = me%nb - 1          !highest index of the basis functions or control points
            allocate(me%kv(0:nk-1), source=k)
                        
            !calculate bezier segments, minimum and maximum regularities
            call me%BSegments(me%nbez,me%unik)
            if(allocated(me%mult)) deallocate(me%mult)
            allocate(me%mult(me%nbez+1))
            do i=1,me%nbez+1
                me%mult(i) = me%FindMult(i=me%FindSpan(me%unik(i)),u=me%unik(i))
            end do
            if(me%nbez>1)then
                me%minreg = me%p - maxval(me%mult(2:me%nbez))
                me%maxreg = me%p - minval(me%mult(2:me%nbez))
            else
                me%minreg = me%p
                me%maxreg = me%p
             end if
            
            !me%ngp = me%get_nquad()

        end subroutine UpdateBasis
        ! .......................................................
            
        !# {!module(basisfunctions)/type(basis)/BSegments.md!}
        pure subroutine BSegments(me,nb,knots)
            class(basis), intent(in)                        :: me
            integer, intent(out)                            :: nb       !< Number of Bezier segments
            real(wp), allocatable, intent(out), optional    :: knots(:) !< The unique knot values
            ! locals
            integer :: i
            real(wp) :: array(me%m-(2*me%p-1))
            
            array = 0.0_wp
            associate(U => me%kv, p => me%p, n => me%n)
                
                ! first, find the number of unique knots
                nb = 1
                array(1) = U(p)
                do i=p+1,n+1
                    if( any( array == U(i) )) cycle
                    nb = nb + 1
                    array(nb) = U(i)
                end do
                
                ! locate the unique knots vector
                if(present(knots))then
                    allocate(knots(nb))
                    knots = array(1:nb)
                end if
                ! substract "1" to determine the number of Bezier segments
                nb = nb - 1
                
            end associate

        end subroutine BSegments
        ! .......................................................
        
        !# {!module(basisfunctions)/type(basis)/ParLen.md!}
        pure function ParLen(me) result(len)
            class(basis), intent(in)    :: me
            real(wp)                    :: len !< Length of the parameter space

            len = me%kv(me%m) - me%kv(0)
        end function ParLen
        ! .......................................................
        
        !# {!module(basisfunctions)/type(basis)/GrevilleAbscissae.md!}
        pure function GrevilleAbscissae(me,dsc) result(pt)
            class(basis), intent(in)        :: me
            logical, optional, intent(in)   :: dsc          !< If true, then the Greville abscissae are modified 
                                                            !< by shifting the boundary points into the domain
            real(wp)                        :: pt(0:me%n)   !< The array containing Greville abscissae
            ! locals
            integer :: i, iak(2), L, li
            real(wp) :: kseq(0:me%n+me%nk), sm(2), ak(2)  

            forall(i=0:me%n) pt(i) = sum(me%kv(i+1:i+me%p))/me%p

            if(present(dsc))then
                if(dsc)then

                kseq(0:me%n)=pt
                kseq(me%n+1:)=me%kv
                kseq = sort(kseq)
                iak(1) = me%p+1; iak(2) = 2*me%nb-1
                L = 1; if(me%p>3) L=2
                sm = 0.0_wp
                do li = 1,L
                    do i=1,2
                        sm(i) = sm(i) + kseq(iak(i)-li) + kseq(iak(i)+li) - 2 * kseq(iak(i))
                    end do
                end do
                ak = sm / (2*L + 1)
                pt(0) = pt(0) + ak(1) 
                pt(me%n) = pt(me%n) + ak(2) 
                end if
            end if

        end function GrevilleAbscissae
        ! .......................................................

        !# {!module(basisfunctions)/type(basis)/get_nquad.md!}
        pure elemental function get_nquad(me,key,reduced_int,pdeq) result(ngp)
            class(basis),intent(in)  :: me
            character(*),intent(in),optional :: key         !< The key defining whether the quadrature is over the Bezier segments or over whole patch  
                                                            !< Valid inputs:  
                                                            !< - ` "elementwise" `     
                                                            !< - ` "patchwise" `  
            logical,intent(in),optional :: reduced_int      !< Reduced integration flag, if true, the routine will return  
                                                            !< the number of quadrature points for reduced integration  
            character(*),intent(in),optional :: pdeq        !< Internal name of the relevant partial differential equation  
                                                            !< Valid inputs:  
                                                            !< - ` "Timoshenko" `  
                                                            !< - ` "Mindlin" `  
                                                            !< - ` "Shell" `  
            !locals
            logical :: r
            integer :: ngp, roi, ngpl

            r=.false.
            
            if(present(pdeq))then
                if((pdeq .iseq. "mindlin") .or. (pdeq .iseq. "shell"))then
                    if(me%p==1)then
                        roi=1; ngpl=2
                    elseif(me%p==2)then
                        roi=2; ngpl=3    
                    elseif(me%p==3)then
                        roi=3; ngpl=3    
                    elseif(me%p==4)then
                        roi=4; ngpl=4      
                    elseif(me%p==5)then
                        roi=4; ngpl=4 
                    elseif(me%p==6)then
                        roi=5; ngpl=5  
                    else
                        ngpl=me%p+1; roi=me%p
                    end if
                elseif(pdeq .iseq. "timoshenko")then
                    if(me%p==1)then
                        roi=1; ngpl=2
                    elseif(me%p==2 .or. me%p==3)then
                        roi=2; ngpl=2   
                    elseif(me%p==4)then
                        roi=3; ngpl=5 
                    elseif(me%p==5 .or. me%p==6 .or. me%p==7 .or. me%p==8)then
                        roi=4; ngpl=5      
                    elseif(me%p==9 .or. me%p==10)then
                        roi=5; ngpl=5
                    else
                        ngpl=me%p+1; roi=me%p
                    end if
                else
                    ngpl=me%p+1; roi=me%p
                end if
            end if

            if(present(reduced_int))then
                r=reduced_int
            end if

            if(present(key))then

                if(key .iseq. "elementwise")then !elementwise

                    ngp=me%p+1
                    if(present(pdeq))then
                        if(r)then
                            ngp = roi
                        else
                            ngp = ngpl
                        end if
                    end if

                elseif(key .iseq. "patchwise")then !patchwise
                    
                    ngp = ceiling(((me%p+1)*me%nbez - (me%minreg+1)*(me%nbez-1))/2.0_wp)
                    if(ngp>64) ngp = 64
                
                end if

            else

                ngp=me%p+1
                if(present(pdeq))then
                    if(r)then
                        ngp = roi
                    else
                        ngp = ngpl
                    end if
                end if

            end if


        end function get_nquad
        ! .......................................................
        
    !***************************************************************************************************************************

    !***************************************************************************************************************************
    ! **TYPE CONSTRUCTORS**  
    ! - basis_constructor   (Fn)  
    ! =>
        
        !# {!module(basisfunctions)/type(basis)/basisConstructor.md!}
        pure function basisConstructor(nk,p,k) result(me)
            integer, intent(in)             :: nk       !< Number of the knots
            integer, intent(in)             :: p        !< Polynomial degree \(p\)
            real(wp), intent(in)            :: k(nk)    !< [Open knot vector](../page/01.fooiga/eqns.html#openvec)  
                                                        !< ` kv(0:m) ` \( = \left\{u_{0},..., u_{m}\right\} \)
            type(basis)                     :: me       !< Spline basis object
            !locals
            integer                     :: i

            me%nk = nk
            me%m  = nk - 1
            me%p  = p
            me%nb  = nk - ( p + 1 )     !number of basis functions or control points
            me%n   = me%nb - 1          !highest index of the basis functions or control points
            allocate(me%kv(0:nk-1), source=k)

            !calculate bezier segments, minimum and maximum regularities
            call me%BSegments(me%nbez,me%unik)
            if(allocated(me%mult)) deallocate(me%mult)
            allocate(me%mult(me%nbez+1))
            do i=1,me%nbez+1
                me%mult(i) = me%FindMult(i=me%FindSpan(me%unik(i)),u=me%unik(i))
            end do
            if(me%nbez>1)then
                me%minreg = me%p - maxval(me%mult(2:me%nbez))
                me%maxreg = me%p - minval(me%mult(2:me%nbez))
            else
                me%minreg = me%p
                me%maxreg = me%p
            end if
                        
        end function basisConstructor
        ! .......................................................


    !***************************************************************************************************************************

    !***************************************************************************************************************************
    ! **PROCEDURES**  
    ! - bin                 (Fn)  
    ! - sort                (Fn)  
    ! =>

        !# {!module(basisfunctions)/type(basis)/bin.md!}
        pure function bin(n, r)
            integer(i8p)        :: bin
            integer, intent(in) :: n
            integer, intent(in) :: r
        
            integer(i8p)        :: num
            integer(i8p)        :: den
            integer             :: i
            integer             :: k
            integer, parameter  :: primes(*) = [2,3,5,7,11,13,17,19]
            num = 1
            den = 1
            do i=0,r-1
                num = num*(n-i)
                den = den*(i+1)
                if (i > 0) then
                    ! Divide out common prime factors
                    do k=1,size(primes)
                        if (mod(i,primes(k)) == 0) then
                            num = num/primes(k)
                            den = den/primes(k)
                        end if
                    end do
                end if
            end do
            bin = num/den
        end function bin
        ! .......................................................
    
        
        !< Recursive quicksort using binary tree pivot.
        pure recursive function sort(x) result(res)
            real(wp), dimension(:), intent(in) :: x !< Input array
            real(wp), dimension(size(x)) :: res
            real(wp), dimension(size(x)-1) :: rest
            real(wp) :: pivot
            if(size(x) > 1)then
                pivot = head(split(x, 2))
                rest = [split(x, 1), tail(split(x, 2))]
                res = [sort(pack(rest, rest < pivot)), pivot, &
                        sort(pack(rest, rest >= pivot))]
            else
                res = x
            endif
        end function sort
    !***************************************************************************************************************************

end module BasisFunctions
