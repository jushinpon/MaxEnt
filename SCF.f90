subroutine SCF (tid,natom,atype,confentropy,atomentropy,filterValue,keepNo,weight,CN_No,CN_ID)
!use information
implicit real*8(a-h,o-z)
integer iscf,ishift 
real*8  oldconfentropy,SCFEntmin,SCFatomentropy(natom)
real*8  oldSCFatomentropy(natom)
real*8  newent,oldent,SCFatkeep(natom)

integer tid,natom,keepeID(natom),CN_No(natom,5),CN_ID(natom,5,50) 
real*8 filterValue,weight(5) 
integer,INTENT(INOUT) :: keepNo,atype(natom)
real*8,INTENT(INOUT) :: confentropy
real*8,INTENT(INOUT) :: atomentropy(natom)

 
!print*,"natom",natom
!print*,"atype",atype
!print*,"confentropy",confentropy
!print*,"atomentropy",atomentropy
!print*,"filterValue",filterValue
!pause

SCFEntmin = confentropy ! set the new value after 
SCFatkeep = atype
SCFatomentropy = atomentropy
newent = 1.		! set different values for the first dowhile
oldent = 1.e20
iscf=0		

do while (oldent .ne. newent)
keepNo = 0
do i = 1,natom
	if(atomentropy(i) .ge. filterValue)then
		keepNo = keepNo + 1
		keepeID(keepNo) = i
	endif	
enddo	 
!print*,"keepeID",keepeID(keepNo)

 iscf = iscf + 1 
 oldent = newent

icheck = 1
   do  ii1=1, keepNo-1           

     do  jj1=ii1+1,keepNo
	 
        ii = keepeID(ii1)
        jj = keepeID(jj1)		
	   
	   if(atype(ii) .ne. atype(jj).and.atomentropy(ii) .ge. filterValue.and.atomentropy(jj).ge.filterValue) then
	      ishift = atype(ii)
          atype(ii) = atype(jj)
          atype(jj) = ishift
		  oldconfentropy = confentropy
		  oldSCFatomentropy = atomentropy ! keep old atom potential
		  
		  temp = confentropy
		
		  call conf_entropyP(ii,jj,natom,atype(:),confentropy,atomentropy(:),weight(:),CN_No(:,:),CN_ID(:,:,:)) !get new confentropy
			if(confentropy .lt. SCFEntmin)then    
			!
			! must work for the first time
				icheck = 0
        SCFEntmin = confentropy
				
				SCFatomentropy = atomentropy				
                SCFatkeep = atype
            else
				ishift = atype(jj)
                atype(jj) = atype(ii)
                atype(ii) = ishift				
                confentropy = oldconfentropy
                atomentropy = oldSCFatomentropy   
				
            endif

       endif
	  
	  end do                      !2 continue
	 
	 end do   	 !1 continue
	 
 atype = SCFatkeep
 atomentropy = SCFatomentropy
 confentropy = SCFEntmin
 newent = SCFEntmin
 
  !if(newent .eq. oldent) then
  !if(icheck .eq. 1 .or. iscf .gt. 30) then !no further decrease by swap
  if(icheck .eq. 1 .or. iscf .gt. 30) then !no further decrease by swap
  
  confentropy = SCFEntmin
  atomentropy = SCFatomentropy
  atype = SCFatkeep

  !write(*,*)""
  exit
  endif 
            
enddo
!print*,"confentropy",confentropy
return
end