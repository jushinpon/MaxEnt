subroutine shuffle
use information
implicit real*8(a-h,o-z)   
integer :: i, randpos, temp
real :: r

do i = size(atype), 2, -1
  call random_number(r)
  randpos = int(r * i) + 1
  temp = atype(randpos)
  atype(randpos) = atype(i)
  atype(i) = temp
end do
  
return 
end
!!!!

recursive subroutine quicksort(a)
  implicit none
  real :: a(:)
  real x, t
  integer :: first = 1, last
  integer i, j

  last = size(a, 1)! one dimension
  x = a( (first+last) / 2 )
  i = first
  j = last
  
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  
  if (first < i - 1) call quicksort(a(first : i - 1))
  if (j + 1 < last)  call quicksort(a(j + 1 : last))
end subroutine quicksort 