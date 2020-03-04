subroutine conf_entropy
use information
implicit real*8(a-h,o-z)   

confentropy = 0.0 ! configurational entropy
atomentropy = 0.0 ! entropy of each atom for normshuffle, set 1 for inital
do 1 i=1,natom
	do 2 ine =1,3 !neighbor atom type
       ntemp = CN_No(i,ine)
		do 3 neID = 1,ntemp ! atom ID of a neighbor type
          JID = CN_ID(i,ine,neID)
           if(ine .le. 2)then !first neighbour atoms
             confentropy = confentropy + weight(ine)*pairweight(atype(i),atype(JID))
             atomentropy(i) = atomentropy(i) + weight(ine)*pairweight(atype(i),atype(JID))
           else
             confentropy = confentropy + weight(ine)
             atomentropy(i) = atomentropy(i) + weight(ine)                     
           endif		  
		  !if(atype(i) .eq. atype(JID))then
          !   confentropy = confentropy + weight(ine)*pairweight(atype(i),atype(JID))
          !   atomentropy(i) = atomentropy(i) + weight(ine)*pairweight(atype(i),atype(JID))
          !endif 
3 		continue                      
2   continue
   ! write(*,*)"atom ",i," atomentropy: ",atomentropy(i)  						
1 continue

!! The following is used for detailed check

!do 123 i=1,natom
!	if(atomentropy(i) .le. 0.0)then
!	
!	write(*,*)"atom: ",i," atomentropy: ",atomentropy(i)
!	
!	do 211 ine =1,2 !neighbor atom type
!       ntemp = CN_No(i,ine)
!		do 311 neID = 1,ntemp ! atom ID of a neighbor type
!          JID = CN_ID(i,ine,neID)
!           if(ine .eq. 1 .or. ine .eq. 2)then !first neighbour atoms
!            ! confentropy = confentropy + weight(ine)*pairweight(atype(i),atype(JID))
!             atomentropy(i) = atomentropy(i) + weight(ine)*pairweight(atype(i),atype(JID))
!			 write(*,*)"neighbourID",ine,"itype",atype(i),"jtype",atype(JID),"jID",JID
!			 write(*,*)"weight",weight(ine),"pairweight",pairweight(atype(i),atype(JID))
!           else
!            ! confentropy = confentropy + weight(ine)
!             atomentropy(i) = atomentropy(i) + weight(ine) 
!             write(*,*)"***neighbourID",ine,"itype",atype(i),"jtype",atype(JID),"jID",JID
!			 write(*,*)"weight",weight(ine)			 
!           endif		  
!		  !if(atype(i) .eq. atype(JID))then
!          !   confentropy = confentropy + weight(ine)*pairweight(atype(i),atype(JID))
!          !   atomentropy(i) = atomentropy(i) + weight(ine)*pairweight(atype(i),atype(JID))
!          !endif 
!311 		continue                      
!211   continue
!
!	write(*,*)''
!	pause
!	endif
!	
!123 continue			
									
!pause
return
end