module Parameters
    
    !< This module declares the global parameters, variables and contains the global procedures

    use, intrinsic :: iso_fortran_env, only : stdout=>output_unit
    
    !the precision parameters and real tolerance values from penf library
    use penf, only: wp => r8p, i8p
    use penf, only: i8p,i4p,i2p,i1p,r16p,r8p,r4p,str
    use penf, only: eps => ZeroR8P 
    use penf, only: smallRPP => smallR8P
    use finer, only: file_ini ! Configuration parsing library
    use vecfor, only: vector

    implicit none
    
    real(wp), parameter      :: pi = 3.1415926535897932384626433832795_wp, pi2 = 6.2831853071795864769252867665590_wp
    real(wp), parameter     :: e_ = 2.7182818284590452353602874713527_wp
    complex(wp), parameter  :: i_ = (0, 1)

    type(file_ini)  :: config     !< Configuration file (*.ini) handler.
    character(256)  :: workspace

    !# Workspace parameters
    character(64)   :: casedef    !< output storage id
    logical         :: clear
    
    !# Settings
    logical         :: quadtreeMethod
    real(wp)        :: accuracy
    logical         :: singularity
    logical         :: plotQuad
    integer         :: nqimg
    
    !# Model parameters
    integer                     :: pd         !< number of parametric dimension
    character(1024)             :: modelfile
    character(64)               :: case
    real(wp)                    :: r1   
    real(wp)                    :: r2   
    type(vector)                :: vpt(3)       !< corners of triangular domain
    integer                     :: npatch
    integer,        allocatable :: ibc(:)
    integer,        allocatable :: bcfaces(:,:)
    character(1),   allocatable :: bctype(:,:)
    integer                     :: ieBCdata
    integer                     :: inBCdata

    !# Convergence parameters
    integer :: p1         !< starting degree
    integer :: npref      !< number of p-ref
    integer :: h0         !< starting value for h-ref
    integer :: nhref      !< number of h-ref
    integer :: dnhref     !< increment of h-ref
    integer :: nm0        !< starting value for mode-ref
    integer :: dnm        !< increment of the number of modes
    integer :: nnm        !< number of mode-ref

    !# Analysis parameters
    character(11)   :: pde   
    character(15)   :: atype      !< analysis type            "staticlinear","eigenanalysis"  
    character(12)   :: itype      !< integration type         !useless     
    character(4)    :: BCond      !< BC type for plates       "1234" S or C, numbers represents sides
    logical         :: reduced
    logical         :: lumped
    integer         :: nmode      !< number of modes

    !# Postprocessing parameters
    character(5)  :: etype      !< external solution type   "exact","FEM"    
    integer       :: res        !< resolution of interpolation

    !# Material parameters
    real(wp) :: El       !< elastic modulus
    real(wp) :: nu       !< poisson ratio
    real(wp) :: Ks       !< shear factor for rectangular sections
    real(wp) :: ro       !< material density
    real(wp) :: Gs       !< shear modulus
    real(wp) :: Cij(5,5)  !< constitutive matrix 

    !# Beam parameters
    !-----------------------------Input----------------------------------------------------!
    real(wp) :: lb       !< beam length
    real(wp) :: wb       !< beam section width
    real(wp) :: hb       !< beam section height
    real(wp) :: ty       !< beam tip load
    real(wp) :: fy       !< beam distributed vertical load
    real(wp) :: mz       !< beam distributed bending moment
    !-----------------------------Derived---------------------------------------------------!
    real(wp) :: Area     !< beam sectional area             != wb*hb
    real(wp) :: Ix       !< beam moment of inertia          != (wb*hb**3)/12
    real(wp) :: EI       !< beam bending rigidity \(EI\)    != El*Ix
    real(wp) :: kGA      !< beam shear rigidity             != Ks*Gs*Area
    
    !# Plate parameters
    !-----------------------------Input----------------------------------------------------!
    real(wp) :: bp       !< plate width
    real(wp) :: lp       !< plate length
    real(wp) :: tp       !< plate thickness
    real(wp) :: fx       !< traction along x-dir (Plane stress)
    ! real(wp) :: fy       !< traction along y-dir (Plane stress)
    real(wp) :: fz       !< plate distributed vertical load
    real(wp) :: mx       !< plate distributed bending moment along x-dir
    real(wp) :: my       !< plate distributed bending moment along y-dir
    !-----------------------------Derived---------------------------------------------------!
    real(wp) :: DSh      !< plate shear rigidity    ! = Ks*tp*Gs                 
    real(wp) :: DB       !< plate flexural rigidty  != El*tp**3/(12*(1-nu**2))   

    !# Plate flow parameters
    real(wp)    :: v0  !< first value of non-dimensional flow velocity
    real(wp)    :: dv  !< increment for non-dimensional flow velocities
    integer     :: nv  !< number of flow velocities
    integer     :: reduce  != 1 !< plate flow ind. vib. eigen key
    real(wp)    :: rhof !< Fluid density
    real(wp)    :: eps2    != 0.001_wp
    logical     :: doublesided
    logical     :: solvedFEM
    logical     :: solvedBEM
    logical     :: integratedHydFrc

    !# Curve or surface fitting option
    integer     :: FittingOpt   !< Keyopt - Interpolation points in parametric coordinates are
                                !<      - 0(default) - the input parameters ` u ` and ` v ` 
                                !<      - 1 - on the control points projection onto patch 
                                !<      - 2 - on the Greville abscissae
                                !<      - 3 - based on the equations given in the NURBS book **Section 9.2.5** 
                                !<      - 4 - based on the equations given in the article by Dornisch et al. (2013)
    integer     :: nip          !< number of interpolation points
    logical     :: check1
    logical     :: check2
    logical     :: check3       !< logical of the checkpoint that visualize the global derivatives

    end module Parameters
    
    !# runs APDL 
    !subroutine apdl_run(file,status)
    !
    !    character(len=*), intent(in) :: file
    !    integer,intent(out) :: status
    !    integer :: lencmd
    !    !character(len=255) fullpath
    !    character(len=32) :: username
    !    character(len=255) :: APDLexe,dir,inpath,outpath
    !    character(len=:), allocatable :: command
    !    character(1024) :: cwd
    !    
    !    call getlog(username)
    !    call getcwd(cwd)
    ! 
    !    ! APDLexe = '"C:\Program Files\ANSYS Inc\v172\ansys\bin\winx64\ansys172.exe"'
    !    APDLexe = '"C:\Program Files\ANSYS Inc\v202\ANSYS\bin\winx64\ANSYS202.exe"'
    !    dir = '"C:\Users\'//trim(username)//'"'
    !    inpath = '"'//trim(cwd)//"\"//trim(file)//'"'
    !    outpath = '"C:\Users\'//trim(username)//'\file.out"'
    !    lencmd = len(APDLexe) + len(dir) + len(inpath) + len(outpath) + 40
    !    allocate (character(len=lencmd) :: command)
    !    command = trim(APDLexe)//' -p ane3fl -dir '//trim(dir)//' -b -i '//trim(inpath)//' -o '//trim(outpath)
    !    call execute_command_line(trim(command),cmdstat=status)
    !    !if (status==0) then
    !    !    write (*,'(a)') "press any key after APDL close itself"
    !    !    pause
    !    !end if
    !    !call runqq(command)
    !    !call system(command)
    !end subroutine apdl_run
    !

    !runtime algorithm
        !integer :: count_0, count_1
        !integer :: count_rate, count_max
        !real(wp) :: time_init, time_final, elapsed_time

        !!Starting time
        !call system_clock(count_0, count_rate, count_max)
        !time_init = count_0*1.0/count_rate
        !do k=1,10000
        !    !do something
        !end do
        !call system_clock(count_1, count_rate, count_max)
        !time_final = count_1*1.0/count_rate
        !! Elapsed time
        !elapsed_time = time_final - time_init
        !write(*,'("  Wall Clock = ",i0,f0.9)') int(elapsed_time),elapsed_time-int(elapsed_time) 
    
    !UNUSED
    
    ! !vecfor precisions
    ! use penf, only : I1P, I2P, I4P, I8P, R4P, R8P, R16P, str
    ! use penf, only : smallRPP=>smallR8P, ZeroRPP=>ZeroR8P !double precision
    !use penf, only : RPP=>R16P, smallRPP=>smallR16P, ZeroRPP=>ZeroR16P !quad precision

    !UIHES precision
    ! use penf, only: wp => R8P, ip => I8P, tiny => smallR8P
    
    ! real(wp), parameter :: thck = 0.01                      !< plate thickness
    ! real(wp), parameter :: DP = elst*thck**3/(12*(1-nu*nu)) !< plate bending rigidity
    ! real(wp), parameter :: kgh   = ks*G*thck                !< 
    ! real(wp), parameter :: C = elst/(1-nu*nu)               !< 
    ! real(wp), parameter :: Cxy = elst/(2-2*nu)              !< 
    
    !output file
    !type(csv_file) :: f
    !logical :: status_ok