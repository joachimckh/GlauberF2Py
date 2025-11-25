module nucmath 
  implicit none
contains

real function pi()
  implicit none
  pi = 4.0 * atan(1.0)
end function pi

subroutine woods_saxon(A, R, d, dmin, coords)
  implicit none
  integer, intent(in) :: A
  real, intent(in) :: R, d, dmin
  real, intent(inout) :: coords(A, 3)
  integer :: i, j
  real :: rr, costh, phi, x, y, z, p, u
  logical :: overlap

  i = 1
  do while (i <= A)
    call random_number(rr)
    rr = 15.0 * rr ! 15 fm max radisu
    call random_number(costh)
    costh = 2.0 * costh - 1.0 ! [-1, 1]
    call random_number(phi)
    phi = 2.0 * pi() * phi
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

subroutine calcPsi2Ecc2(nuc1, nuc2, wounded1, wounded2, A, eps2, psi2)
  implicit none
  integer, intent(in) :: A
  real, intent(in) :: nuc1(A,3), nuc2(A,3)
  real, intent(in) :: wounded1(A), wounded2(A)
  real, intent(out) :: eps2, psi2
  real, allocatable :: P(:,:)
  real, allocatable :: x(:), y(:), r(:), phi(:), w(:)
  real :: cos2, sin2, r2, meanx, meany
  integer :: i, n

  n = sum(wounded1) + sum(wounded2)
  if (n == 0) then
    eps2 = -998.0 
    psi2 = -998.0 
    return
  end if

  allocate(P(n,2))
  P = 0.0

  n = 0
  do i = 1, A
    if (wounded1(i) > 0) then
      n = n + 1
      P(n,1) = nuc1(i,1)
      P(n,2) = nuc1(i,2)
    end if
  end do
  do i = 1, A
    if (wounded2(i) > 0) then
      n = n + 1
      P(n,1) = nuc2(i,1)
      P(n,2) = nuc2(i,2)
    end if
  end do

  allocate(x(n), y(n), r(n), phi(n), w(n))

  meanx = sum(P(:,1)) / real(n)
  meany = sum(P(:,2)) / real(n)

  x = P(:,1) - meanx
  y = P(:,2) - meany

  r = sqrt(x**2 + y**2)
  phi = atan2(y, x)
  w = r**2

  cos2 = sum(w * cos(2.0 * phi))
  sin2 = sum(w * sin(2.0 * phi))
  r2 = sum(w)

  if (r2 == 0.0) then
    eps2 = -999.0 
    psi2 = -999.0 
    deallocate(P, x, y, r, phi, w)
    return
  end if

  eps2 = sqrt(cos2**2 + sin2**2) / r2
  psi2 = 1/2.0 * (atan2(sin2, cos2) + pi())

  deallocate(P, x, y, r, phi, w)
end subroutine calcPsi2Ecc2



end module nucmath 






