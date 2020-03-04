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
! making different configuration by shuffle or other ways

! the following is for the case without any better case found
! after the do while 

SCFatkeep = atype
SCFatomentropy = atomentropy



newent = 1.		! set different values for the first dowhile
oldent = 1.e20
iscf=0		
!write(*,*)"SCF step: ", 0," current entropy: ",SCFEntmin/dble(natom)


do while (oldent .ne. newent)
!!!!!!!!!!
keepNo = 0
do i = 1,natom
	if(atomentropy(i) .ge. filterValue)then
		keepNo = keepNo + 1
		
		keepeID(keepNo) = i
		!write(*,*)i,atomentropy(i),keepNo
		!pause
	endif	
	!if(keepNo .ge. 200) exit ! only consider the maximal number of 100
enddo
!write(*,*)i,atomentropy(i),keepNo     
	 
if(keepNo .eq. 0)then
 write(*,*)"***All atoms with larger scoring values have been processed!!***"
 write(*,*)"YOU CAN TERMINATE YOUR CODE or WAIT for BETTER by MC!!!!!!!!!!!"
 
 !do i =1,natom
 !   write(*,*)i,atomentropy(i),keepNo
 !enddo
 !  write(*,*)"Ave confentropy:", confentropy/dble(natom) 
 !pause
endif 
!!!!!!!!!!!!!!!!!!!

 iscf = iscf + 1 
 oldent = newent

  

!$OMP PARALLEL &
!$OMP shared(SCFEntmin,SCFatkeep,SCFatomentropy,atype) &
!$OMP private(jj,ishift,oldSCFatomentropy,oldconfentropy)&
!$OMP firstprivate(confentropy,atomentropy,natom)

!write(*,*)"SCF time",iscf,"Keepno",keepNo
!write(*,*)"Before thread ID:",omp_get_thread_num()
!!write(*,*)"SCFatomentropy:",SCFatomentropy
!write(*,*)"Before SCFEntmin",SCFEntmin
!write(*,*)"=="

!write(*,*)iscf,"Before thread ID:",omp_get_thread_num(),SCFEntmin

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
		  !write(*,*)"****oldconfentropy",oldconfentropy,ii,jj
		  !write(*,*)"T", omp_get_thread_num(),ii,jj
		  !write(*,*)confentropy
		  temp = confentropy
		  call conf_entropyP(ii,jj) !get new confentropy
		  !call conf_entropy
		 ! if(confentropy .lt.0)then
		 !  write(*,*)"T", omp_get_thread_num(),ii,jj
		 !  write(*,*)confentropy,temp
		 !  
		 !  pause
		 !  exit
		 ! 
		 ! endif
		  !write(*,*)"****afterP",confentropy
		  !pause
!
			if(confentropy .lt. SCFEntmin)then     !小於
			!
			! must work for the first time
			    !$OMP CRITICAL 
				icheck = 0
                SCFEntmin = confentropy
				
				SCFatomentropy = atomentropy				
                !xmin = x
                !ymin = y
                !zmin = z
				
                SCFatkeep = atype
		       !$OMP END CRITICAL 
            else
				
                !swap back
				ishift = atype(jj)
                atype(jj) = atype(ii)
                atype(ii) = ishift				
				!x = xmin
                !y = ymin
                !z = zmin
				!atype = atkeep
                confentropy = oldconfentropy
                atomentropy = oldSCFatomentropy   
				
            endif
!!			
       endif
	  
	  end do                      !2 continue
	 
	 end do   	 !1 continue
	 
!write(*,*)"keepNo",keepNo,"SCF",iscf, "Ave Ent",confentropy/dble(natom)
!$OMP END do
!!$omp barrier
!write(*,*)"SCF time",iscf
!write(*,*)iscf,"thread ID:",omp_get_thread_num(),SCFEntmin
!write(*,*)"SCFatomentropy:",SCFatomentropy
!write(*,*)"SCFEntmin",SCFEntmin
!write(*,*)""
!write(*,*)"**************"

!$OMP END PARALLEL

! x = xmin
! y = ymin
! z = zmin
!! use the best parameters for the next while loop
!write(*,*)"SCF step: ", iscf," current entropy: ",SCFEntmin/dble(natom)

! if no lower value is found after the loop, it is ok because the initial values for SCF parameters
! were assigned by the original related values.
 atype = SCFatkeep
 atomentropy = SCFatomentropy
 confentropy = SCFEntmin
 newent = SCFEntmin
 
  !if(newent .eq. oldent) then
  if(icheck .eq. 1 .or. iscf .gt. 30) then !no further decrease by swap
  !write(*,*)""
  !write(*,*)"*****current best  entropy after SCF process",SCFEntmin/dble(natom)
  !write(*,*)"####CURRENT system best entropy", Entmin/dble(natom)
  !write(*,*)"Exit, no more better results!!"
  !put SCF results to the corresponding global parameters
  confentropy = SCFEntmin
  atomentropy = SCFatomentropy
  atype = SCFatkeep

  !write(*,*)""
  exit
  endif 
            
enddo


 

return
end