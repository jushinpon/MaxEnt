subroutine twoshuffle
use information
implicit real*8(a-h,o-z)   
integer i1,i2,ishift
real*8 :: r,temp,tempent

! pick the first one
	call random_number(r)
    i1 = int(natom*r) + 1
! pick the second one  
1112 continue
	call random_number(r)
    i2 = int(natom*r) + 1
	if(i1 .eq.i2) goto 1112
! start the shift
    ishift = atype(i1)
    atype(i1) = atype(i2)
    atype(i2) = ishift     

return
end