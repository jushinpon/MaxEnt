subroutine conf_entropyP(ii,jj)
use information
implicit real*8(a-h,o-z)   
integer ii,jj ! the swapped atoms 
integer::ikeep(500),ipair(2)

ikeep = 0
ipair(1) = ii
ipair(2) = jj
!print*, "In con_entropyP ii jj: ",ii,jj
!		  print *, "confentropy not updated type switched: ",confentropy/dble(natom)
		  
!write(*,*)"****+++ enter sub confentropy",confentropy
!confentropy = 0.0 ! configurational entropy
!atomentropy = 1.0 ! entropy of each atom for normshuffle, set 1 for inital

ndup = 0 ! atom ID counter without duplicate ID

do 1 ip=1,2
	i = ipair(ip)
!check two switched atoms first
		if(ndup .eq.0)then
			ndup = ndup + 1
			ikeep(ndup) = i ! keep the id for new atomentropy
		else
			index = 0
				do idu = 1,ndup ! check if this atom ID has been kept
						if(ikeep(idu) .eq. i)then
						index = 1
						exit
						endif
				enddo
				if(index .eq. 0)then !no duplicate atom ID
						ndup = ndup + 1
					if(ndup .gt. 500) write(*,*)"!!ikeep array dead!!"
						ikeep(ndup) = i
					endif			 
				endif

!!the following is to analyze the neighbour atoms of two ipair atoms
	do 2 ine =1,1 !neighbor atom type
       ntemp = CN_No(i,ine)
		do 3 neID = 1,ntemp ! atom ID of a neighbor type
          JID = CN_ID(i,ine,neID) ! the nth neighbor atom ID of i
		  
		  !if(ndup .eq.0)then
		  !  ndup = ndup + 1
			!ikeep(ndup) = JID ! keep the id for new atomentropy
		  !else
		      index = 0
		     do idu = 1,ndup ! check if this atom ID has been kept
		        if(ikeep(idu) .eq. JID)then
				  index = 1
				  exit
				endif
		     enddo
		       if(index .eq. 0)then !no duplicate atom ID
			      ndup = ndup + 1
				    if(ndup .gt. 500) write(*,*)"!!ikeep array dead!!"
			      ikeep(ndup) = JID
			   endif			 
		  !endif		 
3 		continue                      
2   continue  						
1 continue

pconfentropy = 0.0 ! need to be deducted from original confentropy

do 111 ip=1,ndup
  i = ikeep(ip)
  pconfentropy = pconfentropy + atomentropy(i) ! get the sum of original partial 
111 continue

!write(*,*)"ii",ii,"jj",jj
!write(*,*)"original confentropy",confentropy

confentropy = confentropy - pconfentropy ! deduct those changed after swap

 
do 11 ip=1,ndup
  i = ikeep(ip)
  atomentropy(i) = 0.0
	do 21 ine =1,1 !neighbor atom type
       ntemp = CN_No(i,ine)
		do 31 neID = 1,ntemp ! atom ID of a neighbor type
          JID = CN_ID(i,ine,neID) 
		!  if(ine .le. 1 .and. atype(i) .eq. atype(JID))then !first neighbour atoms
        !     atomentropy(i) = atomentropy(i) + weight(ine)*pairweight(atype(i),atype(JID))
        !     
        !   else
             if(atype(i) .eq. atype(JID))then
				atomentropy(i) = atomentropy(i) + weight(ine)
             endif                     
        !   endif
31 		continue                      
21   continue
  	confentropy = confentropy + atomentropy(i)					
11 continue

!write(*,*)"ii jj: ",ii,jj
!write(*,*)"***Before conf_entropy",confentropy/dble(natom)
!call conf_entropy
!write(*,*) "***After conf_entropy: ",confentropy/dble(natom)
!print*, "switch back"
!ishift = atype(ii)
!atype(ii) = atype(jj)
!atype(jj) = ishift
!call conf_entropy
!write(*,*) "***After switch: ",confentropy/dble(natom)
!	
!pause

return
end
