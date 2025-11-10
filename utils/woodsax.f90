module sampler
  implicit none
contains

subroutine woods_saxon(A, R, d, dmin, coords)
  implicit none
  integer, intent(in) :: A
  real, intent(in) :: R, d, dmin
  real, intent(inout) :: coords(A, 3)
  real, parameter :: pi = 3.141592653589793
  integer :: i, j
  real :: rr, costh, phi, x, y, z, p, u
  logical :: overlap

  i = 1
  do while (i <= A)
    call random_number(rr)
    rr = 15.0 * rr
    call random_number(costh)
    costh = 2.0 * costh - 1.0
    call random_number(phi)
    phi = 2.0 * pi * phi
    p = 1.0 / (1.0 + exp((rr - R) / d))
    call random_number(u)
    if (u < p) then
      x = rr * sqrt(1.0 - costh**2) * cos(phi)
      y = rr * sqrt(1.0 - costh**2) * sin(phi)
      z = rr * costh
      overlap = .false.
      do j = 1, i - 1
        if (sqrt((x - coords(j,1))**2 + (y - coords(j,2))**2 + (z - coords(j,3))**2) < 2.0 * dmin) then
          overlap = .true.
          exit
        end if
      end do
      if (.not. overlap) then
        coords(i,1) = x
        coords(i,2) = y
        coords(i,3) = z
        i = i + 1
      end if
    end if
  end do
end subroutine woods_saxon

end module sampler
