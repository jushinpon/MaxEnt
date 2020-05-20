subroutine kshuffle(tid,imc,natom,atype,confentropy,atomentropy,weight,CN_No,CN_ID)
!use information
implicit real*8(a-h,o-z)   
integer i1,i2,ishift,imc,iposNo
real*8 r,temp,tempent,tempconfentropy
integer positiveID(natom) 

integer tid,nkMCshift,natom
real*8 weight(5)
integer CN_No(natom,5),CN_ID(natom,5,50)
integer,INTENT(INOUT) :: atype(natom)
real*8,INTENT(INOUT) :: confentropy
real*8,INTENT(INOUT) :: atomentropy(natom)

call random_number(r)
!nkMCshift gradually becomes smaller with the increasing MC steps
!if(mod(imc,20).eq.0)then
nkMCshift = dint(natom*0.05*r)+1
!nkMCshift = 10
if(nkMCshift .lt. 10) nkMCshift = 10
!write(*,*)"nkMCshift 1:",nkMCshift,imc
!pause
!else
!nkMCshift = int(20*r)+1 ! the times for kMC-like shift, maximum 20 times
!write(*,*)"nkMCshift 2:",nkMCshift,imc
!endif

!nkMCshift = 1
!print*, "in kshuffle before do 1, nkMCshift",nkMCshift
do 1 ikMC=1, nkMCshift

	call random_number(r)	
	
    temp = r*confentropy
    tempent = 0.0
 !   print*, "in kshuffle before natom do 1 iposNo: ",iposNo
    i1 = 0
    	do i = 1,natom
        	!print*, "in natom loop,ip: ",ip

			!i = positiveID(ip)
        	!print*, "in natom loop,i: ",i

        	tempent = tempent +  atomentropy(i)
        	!print*, "in natom loop,i and atomentropy(i): ",i,atomentropy(i)
        	!print*, "in natom loop: temp,tempent: ",temp,tempent

            if(tempent .gt.temp)then
               i1 = i
           	!print*, "****in natom IF loop: i and i1: ",i,i1

               exit
            endif   
        enddo
if (i1 .gt. 0)then  ! i1 has been assiged a ID from above      
!print*, "in kshuffle after natom do 2 i1: ",i1
irepeat  = 0 ! parameter to count the 1112 continuoue times
! pick the second one  
1112 continue

irepeat  = irepeat + 1
if (irepeat .gt. 200) goto 12345 ! go to randomly pick another atoms 
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
    !i2 = positiveID(i2temp)
	!.or. atomentropy(i2) .ge. filterValue
!	    print*, "in kshuffle before goto 1112 1",i1,i2

	if(i1 .eq. i2 .or. atype(i1) .eq. atype(i2)) goto 1112
! start the shift
!	    print*, "in kshuffle before after 1112 1",i1,i2

    ishift = atype(i1)
    atype(i1) = atype(i2)
    atype(i2) = ishift  
 !   print*, "in kshuffle before conf_entropyP",i1,i2

  call conf_entropyP(i1,i2,natom,atype,confentropy,atomentropy,weight,CN_No,CN_ID)     
  ! print*, "in kshuffle After_entropyP",i1,i2


else ! all atomentropy(i) are negative or goto 1112 over 100 times

12345 continue ! if goto 1112 over 200 times not work 

	call random_number(r)
  
	i1pos = int(r * natom) + 1
11111 continue	
	call random_number(r)
  
	i2pos = int(r * natom) + 1
	
	if(i1pos .eq. i2pos .or. atype(i1pos) .eq. atype(i2pos)) goto 11111

	itemp = atype(i1pos)
	atype(i1pos) = atype(i2pos)
	atype(i2pos) = itemp
    call conf_entropyP(i1pos,i2pos,natom,atype,confentropy,atomentropy,weight,CN_No,CN_ID)
endif

1 continue ! atom type swap loop

!print*, "in kshuffle after do 2"

return
end
