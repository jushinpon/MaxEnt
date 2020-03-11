subroutine kshuffle(imc)
use information
implicit real*8(a-h,o-z)   
integer i1,i2,ishift,imc
real*8 :: r,temp,tempent

call random_number(r)
!nkMCshift gradually becomes smaller with the increasing MC steps
!if(mod(imc,20).eq.0)then
nkMCshift = int(natom*0.02*r)
!nkMCshift = 10
if(nkMCshift .lt. 2) nkMCshift = 6
!write(*,*)"nkMCshift 1:",nkMCshift,imc
!pause
!else
!nkMCshift = int(20*r)+1 ! the times for kMC-like shift, maximum 20 times
!write(*,*)"nkMCshift 2:",nkMCshift,imc
!endif

!nkMCshift = 1
do 1 ikMC=1, nkMCshift
! pick the first one  
	call random_number(r)
    temp = r*confentropy
    tempent = 0.0
    	do i = 1,natom
        	tempent = tempent +  atomentropy(i)
            if(tempent .gt.temp)then
               i1 = i
               exit
            endif   
        enddo

! pick the second one  
1112 continue
	call random_number(r)
!    temp = r*confentropy
!    tempent = 0.0
!    	do i = 1,natom
!        	tempent = tempent +  atomentropy(i)
!            if(tempent .gt.temp)then
!               i2 = i
!               exit
!            endif   
!        enddo
    i2 = int(natom*r)+1
	!.or. atomentropy(i2) .ge. filterValue
	if(i1 .eq. i2 .or. atype(i1) .eq. atype(i2)) goto 1112
! start the shift
    ishift = atype(i1)
    atype(i1) = atype(i2)
    atype(i2) = ishift  
  call conf_entropyP(i1,i2)     
1 continue

return
end
