subroutine randshuffle
use information
implicit real*8(a-h,o-z)   
integer i1,i2,ishift
real*8 :: r,temp,tempent

call random_number(r)
!nkMCshift gradually becomes smaller with the increasing MC steps
nkMCshift = int(natom*acceptratio*r)+1 ! the times for kMC-like shift
nkMCshift = 1
do 1 ikMC=1, nkMCshift
! pick the first one
	call random_number(r)
    i1 = int(natom*r) + 1
! pick the second one  
1112 continue
	call random_number(r)
    i2 = int(natom*r) + 1
	if(i1 .eq.i2 .or. atype(i1) .eq. atype(i2)) goto 1112
! start the shift
    ishift = atype(i1)
    atype(i1) = atype(i2)
    atype(i2) = ishift     
1 continue

return
end