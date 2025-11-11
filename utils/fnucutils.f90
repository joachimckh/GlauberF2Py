
module nucutils
  use nucmath 
  implicit none

contains



subroutine find_collisions(nucleiA, nucleiB, sig_nn, size, part1, part2)
  implicit none

  integer, intent(in) :: size
  real, intent(in) :: nucleiA(size,3), nucleiB(size,3)
  real, intent(in) :: sig_nn
  ! logical, intent(inout) :: part1(size), part2(size) ! doesnt work with f2py..
  real, intent(inout) :: part1(size), part2(size)
  ! real, intent(inout) :: bc(size,3) ! not needed for the purposes currently

  integer :: i, j
  real :: dx, dy, dz, dist 
  real :: sig_nn_eff

  sig_nn_eff = sig_nn / pi() 

  ! part1 = .false.
  ! part2 = .false.

  do i = 1, size
    do j = 1, size
      dx = nucleiA(i,1) - nucleiB(j,1)
      dy = nucleiA(i,2) - nucleiB(j,2)
      dz = nucleiA(i,3) - nucleiB(j,3)
      dist = dx*dx + dy*dy + dz*dz
      if (dist < sig_nn_eff) then
        part1(i) = 1!.true. 
        part2(j) = 1!.true. 
        ! bc(i,1) = (nucleiA(i,1) + nucleiB(j,1) )/ 2.0
        ! bc(i,2) = (nucleiA(i,2) + nucleiB(j,2) )/ 2.0
        ! bc(i,3) = (nucleiA(i,3) + nucleiB(j,3) )/ 2.0
      end if
    end do
  end do

end subroutine find_collisions

end module nucutils 






