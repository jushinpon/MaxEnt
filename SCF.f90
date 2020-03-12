subroutine SCF

use information
implicit real*8(a-h,o-z)
integer ( kind = 4 ) tid
integer iscf,ishift 
real*8  oldconfentropy,SCFEntmin,SCFatomentropy(natom)
real*8  oldSCFatomentropy(natom)
real*8  newent,oldent,SCFatkeep(natom)

!pause
!!!!!call omp_set_num_threads(1)

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

 iscf = iscf + 1 
 oldent = newent

!$OMP PARALLEL &
!$OMP shared(SCFEntmin,SCFatkeep,SCFatomentropy,atype) &
!$OMP private(jj,ishift,oldSCFatomentropy,oldconfentropy)&
!$OMP firstprivate(confentropy,atomentropy,natom)

!$OMP do
icheck = 1
   do  ii1=1, keepNo-1           

     do  jj1=ii1+1,keepNo
	 
        ii = keepeID(ii1)
        jj = keepeID(jj1)		
	   
	   if(atype(ii) .ne. atype(jj).and.atomentropy(ii) .ge. filterValue.and.atomentropy(jj).ge.filterValue) then	        !ne 不等於
	      ishift = atype(ii)
          atype(ii) = atype(jj)
          atype(jj) = ishift
		  oldconfentropy = confentropy
		  oldSCFatomentropy = atomentropy ! keep old atom potential
		  
		  temp = confentropy
		
		  call conf_entropyP(ii,jj) !get new confentropy
			if(confentropy .lt. SCFEntmin)then     !小於
			!
			! must work for the first time
			    !$OMP CRITICAL 
				icheck = 0
                SCFEntmin = confentropy
				
				SCFatomentropy = atomentropy				
                SCFatkeep = atype
		       !$OMP END CRITICAL 
            else
				ishift = atype(jj)
                atype(jj) = atype(ii)
                atype(ii) = ishift				
                confentropy = oldconfentropy
                atomentropy = oldSCFatomentropy   
				
            endif
!!			
       endif
	  
	  end do                      !2 continue
	 
	 end do   	 !1 continue
	 
!$OMP END do

!$OMP END PARALLEL

 atype = SCFatkeep
 atomentropy = SCFatomentropy
 confentropy = SCFEntmin
 newent = SCFEntmin
 
  !if(newent .eq. oldent) then
  if(icheck .eq. 1 .or. iscf .gt. 30) then !no further decrease by swap
  
  confentropy = SCFEntmin
  atomentropy = SCFatomentropy
  atype = SCFatkeep

  !write(*,*)""
  exit
  endif 
            
enddo

return
end
