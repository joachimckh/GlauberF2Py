
module simulation
  use nucmath
  use nucutils
  implicit none

contains

  subroutine proc(n_events, b_array, n_b, A, results, psi2, ecc2)
    implicit none
    integer, intent(in) :: n_events
    real, intent(in) :: b_array(n_b)
    integer, intent(in) :: n_b
    integer, intent(in) :: A 

    real, intent(inout) :: results(n_events,n_b)
    real, intent(inout) :: psi2(n_events, n_b), ecc2(n_events, n_b)
    
    integer :: i_event, i_b
    ! real, allocatable :: ecc2_array(:), psi2_array(:)
    
    real :: R, d, dmin
    real :: nuc1( A, 3 )
    real :: nuc2( A, 3 )
    real :: part1( A )
    real :: part2( A )

    real :: bshift
    real :: eps2, vpsi2

    !allocate(ecc2_array(n_b))
    !allocate(psi2_array(n_b))

    

    R = 1.25 * A**(1.0/3.0)
    d = 0.5
    dmin = 0.4

    do i_b = 1, n_b
      bshift = b_array(i_b) / 2.0
      do i_event = 1, n_events
        call woods_saxon(A,R,d, dmin, nuc1)
        call woods_saxon(A,R,d, dmin, nuc2)
        nuc1(:,1) = nuc1(:,1) - bshift
        nuc2(:,1) = nuc2(:,1) + bshift
        call find_collisions(nuc1, nuc2, 10.0, A, part1, part2)
        call calcPsi2Ecc2(nuc1, nuc2, part1, part2, A, eps2, vpsi2) 
        ! count collisions
        results(i_event, i_b) = sum(part1) + sum(part2) 
        ecc2(i_event, i_b ) = eps2 
        psi2(i_event, i_b) = vpsi2
        part1 = 0.0
        part2 = 0.0
      end do
    end do




    ! deallocate(ecc2_array)
    ! deallocate(psi2_array)

  end subroutine proc 


end module simulation
