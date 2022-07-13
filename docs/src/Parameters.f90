module Parameters
    
    !< This module declares the global parameters, variables and contains the global procedures

    use, intrinsic :: iso_fortran_env, only : stdout=>output_unit
    
    !the precision parameters and real tolerance values from penf library
    use penf, only: wp => r8p, i8p
    use penf, only: eps => ZeroR8P 

    implicit none
    
    real(wp), parameter      :: pi = 3.1415926535897932384626433832795_wp, pi2 = 6.2831853071795864769252867665590_wp
    real(wp), parameter     :: e_ = 2.7182818284590452353602874713527_wp
    complex(wp), parameter  :: i_ = (0, 1)


    end module Parameters
    