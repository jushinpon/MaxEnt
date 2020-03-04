subroutine partialshuffle
use information
implicit real*8(a-h,o-z)   
integer :: i, randpos, temp,itemp1,itemp2,i1,i2
real :: r

!npart = 5 ! maximal atoms to change types
!itemp = natom - npart
!
!111 continue
!
!call random_number(r)
!itemp1 = int(r * dble(natom) ) + 1
! 
!if(itemp1 .gt. npart .and. itemp1 .lt. itemp)then
!    call random_number(r)
!	itemp2 = itemp1 + int(2.*dble(npart)*r - npart) ! -npart to + npart
!!	write(*,*)'case 1'
!elseif(itemp1 .le. npart)then
!    itemp2 = int(dble(itemp1-1)*r)+1
!!	write(*,*)'case 2'
!elseif(itemp1 .ge. itemp)then  
!    itemp2 = itemp1 + int(dble(natom-itemp1)*r)+1
!!	write(*,*)'case 3'
!endif

111 continue

call random_number(r)
itemp1 = int(r * dble(natom) ) + 1
call random_number(r)
itemp2 = int(r * natom) + 1

if(itemp1 .eq. itemp2) then
  goto 111
endif

!
if(itemp1 .gt. itemp2)then
  i1 =  itemp1
  i2 =  itemp2
elseif(itemp1 .lt. itemp2) then
  i1 =  itemp2
  i2 =  itemp1
elseif(itemp1 .eq. itemp2) then
  goto 111
endif

!temp = atype(itemp1)
!atype(itemp1) = atype(itemp2)
!atype(itemp2) = temp

do i = i1, i2, -1
  call random_number(r)
  randpos = int(natom * r) + 1
  temp = atype(randpos)
  atype(randpos) = atype(i)
  atype(i) = temp
end do
  
return 
end 