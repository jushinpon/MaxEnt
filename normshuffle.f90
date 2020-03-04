subroutine normshuffle
use information
implicit real*8(a-h,o-z)   
integer i1,i2,ishift
real*8 :: r,temp,tempent,inverEntsum,inveratomentropy(natom)

inverEntsum = 0.0
inveratomentropy = 0.0

do i=1,natom
    if(atomentropy(i) .eq. 0) atomentropy(i) = 1e-3
	inveratomentropy(i) = 1.d0/atomentropy(i)
	inverEntsum = inverEntsum + inveratomentropy(i)
enddo

call random_number(r)
!nkMCshift gradually becomes smaller with the increasing MC steps
!nkMCshift = int(natom*acceptratio*r)+1 ! the times for kMC-like shift
nkMCshift = 1
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

! pick the second one  (use the inverse of entropy)
nsecond = 0
1112 continue
	call random_number(r)
	nsecond = nsecond +1
    temp = r*inverEntsum
    tempent = 0.0
    	do i = 1,natom
        	tempent = tempent +  inveratomentropy(i)
            if(tempent .gt.temp)then
               i2 = i
			   exit
            endif   
        enddo
		
	if(i1 .eq. i2 .or. atype(i1) .eq. atype(i2)) then
	if(inverEntsum .gt.1000) then
	write(*,*)" time to find second atom: ",nsecond
	write(*,*)"random number: ",r, temp
	write(*,*)"i1= ",i1
	write(*,*)atomentropy(i1),inveratomentropy(i1)
	write(*,*)"i2= ",i2
	write(*,*)atomentropy(i2),inveratomentropy(i2)
	
	write(*,*)"type i1= ",atype(i1)
	write(*,*)"type i2= ",atype(i2)
	write(*,*)"inverEntsum: ",inverEntsum
	write(*,*)"tempent ",tempent
	endif
!if(inverEntsum .gt.1000) then
! do inv = 1,natom
!   write(*,*)"inv: ",inv, inveratomentropy(inv)
! enddo
!endif
!pause	
	goto 1112
	endif
! start the shift
    ishift = atype(i1)
    atype(i1) = atype(i2)
    atype(i2) = ishift  
  !call conf_entropy     
1 continue

return
end