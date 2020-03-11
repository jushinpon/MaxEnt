subroutine conf_entropyP(ii,jj)
use information
implicit real*8(a-h,o-z)   
integer ii,jj ! the swapped atoms 
integer::ikeep(500),ipair(2)

ikeep = 0
ipair(1) = ii
ipair(2) = jj

!write(*,*)"****+++ enter sub confentropy",confentropy
!confentropy = 0.0 ! configurational entropy
!atomentropy = 1.0 ! entropy of each atom for normshuffle, set 1 for inital

ndup = 0 ! atom ID counter without duplicate ID

do 1 ip=1,2
  i = ipair(ip)
	do 2 ine =1,3 !neighbor atom type
       ntemp = CN_No(i,ine)
		do 3 neID = 1,ntemp ! atom ID of a neighbor type
          JID = CN_ID(i,ine,neID) ! the nth neighbor atom ID of i
		  if(ndup .eq.0)then
		    ndup = ndup + 1
			ikeep(ndup) = JID ! keep the id for new atomentropy
		  else
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
		  endif		 
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
	do 21 ine =1,3 !neighbor atom type
       ntemp = CN_No(i,ine)
		do 31 neID = 1,ntemp ! atom ID of a neighbor type
          JID = CN_ID(i,ine,neID) 
		  if(ine .le. 1 )then !first neighbour atoms
            ! confentropy = confentropy + weight(ine)*pairweight(atype(i),atype(JID))
             atomentropy(i) = atomentropy(i) + weight(ine)*pairweight(atype(i),atype(JID))
           else
             !confentropy = confentropy + weight(ine)
             atomentropy(i) = atomentropy(i) + weight(ine)                     
           endif
		 ! if(atype(i) .eq. atype(JID))then             
         !    atomentropy(i) = atomentropy(i) + weight(ine)*pairweight(atype(i),atype(JID))
         ! endif 
31 		continue                      
21   continue
  	confentropy = confentropy + atomentropy(i)					
11 continue

!write(*,*)"pconfentropy",pconfentropy,"modified",confentropy
!write(*,*)"ndp",ndup	
	
!pause

return
end
