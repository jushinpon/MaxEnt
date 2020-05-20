subroutine shuffle(natom,atype)
implicit real*8(a-h,o-z)   
!use information
integer natom
integer atype(natom)
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