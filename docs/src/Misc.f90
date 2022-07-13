module misc
    
    use Parameters, only:wp, eps
    
    implicit none

    interface operator(.iseq.)
        module procedure real_equal
        module procedure char_equal
    end interface

    contains
    
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
    
end module misc