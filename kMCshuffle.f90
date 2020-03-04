subroutine kMCshuffle
use information
implicit real*8(a-h,o-z)   
integer i1,i2,ishift
real*8 :: r,temp,tempent

call random_number(r)
ndivid = int(r*10) + 2 ! min = 2, max = 11
nkMCshift = int(natom/ndivid) ! the times for kMC-like shift

do 1 ikMC=1, nkMCshift
! pick the first one  
	call random_number(r)
    temp = r*confentropy
    tempent = 0.0
    	do i = 1,natom
        	tempent = tempent +  atomentropy(i)
            if(tempent .gt.temp .and. i.eq.1) then
              i1 =1
              exit
            elseif(tempent .gt.temp .and. i.ne.1)then
              i1 = i-1
              exit
            endif  
        enddo

! pick the second one  
1112 continue
	call random_number(r)
    temp = r*confentropy
    tempent = 0.0
    	do i = 1,natom
        	tempent = tempent +  atomentropy(i)
            if(tempent .gt.temp .and. i.eq.1)then
               i2 =1
               exit
            elseif(tempent .gt.temp .and. i.ne.1)then
               i2 = i-1
               exit
            endif   
        enddo
	if(i1 .eq. i2) goto 1112
! start the shift
    ishift = atype(i1)
    atype(i1) = atype(i2)
    atype(i2) = ishift     
1 continue
return
end