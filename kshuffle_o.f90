subroutine kshuffle(imc)
use information
implicit real*8(a-h,o-z)   
integer i1,i2,ishift,imc
real*8 :: r,temp,tempent

call random_number(r)
!nkMCshift gradually becomes smaller with the increasing MC steps
!if(mod(imc,20).eq.0)then
nkMCshift = dint(natom*0.02*r)+1
!nkMCshift = 10
if(nkMCshift .le. 2) nkMCshift = 6
!write(*,*)"nkMCshift 1:",nkMCshift,imc
!pause
!else
!nkMCshift = int(20*r)+1 ! the times for kMC-like shift, maximum 20 times
!write(*,*)"nkMCshift 2:",nkMCshift,imc
!endif

!nkMCshift = 1
print*, "in kshuffle before do 1, nkMCshift",nkMCshift
do 1 ikMC=1, nkMCshift
! pick the first one  
	call random_number(r)
	
	
	
    temp = r*confentropy
    tempent = 0.0
    print*, "in kshuffle before natom do 1"

    	do i = 1,natom
        	tempent = tempent +  atomentropy(i)
        	print*, "in natom loop,i and atomentropy(i): ",i,atomentropy(i)
        	print*, "in natom loop: temp,tempent: ",temp,tempent

            if(tempent .gt.temp)then
               i1 = i
           	print*, "****in natom IF loop: i and i1: ",i,i1

               exit
            endif   
        enddo
print*, "in kshuffle after natom do 2 i1: ",i1

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
    i2 = dint(natom*r)+1
	!.or. atomentropy(i2) .ge. filterValue
	    print*, "in kshuffle before goto 1112 1",i1,i2

	if(i1 .eq. i2 .or. atype(i1) .eq. atype(i2)) goto 1112
! start the shift
	    print*, "in kshuffle before after 1112 1",i1,i2

    ishift = atype(i1)
    atype(i1) = atype(i2)
    atype(i2) = ishift  
    print*, "in kshuffle before conf_entropyP",i1,i2

  call conf_entropyP(i1,i2)     
   print*, "in kshuffle After_entropyP",i1,i2

1 continue
print*, "in kshuffle after do 2"

return
end
