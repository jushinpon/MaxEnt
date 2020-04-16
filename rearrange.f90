subroutine rearrange
use information
implicit real*8(a-h,o-z)   
!for rearrange types
integer,allocatable::temptype(:),occupied(:),CN_typeNo(:,:),minimumIDs(:)
integer,allocatable::unoccupied(:)
integer unoccupiedCount,itempid(2),minimumCounter
!integer::temptype(natom),occupied(natom),CN_typeNo(natom,elemtype),minimumIDs(natom)
!integer::unoccupied(natom)
!integer unoccupiedCount,itempid(2),minimumCounter

allocate(temptype(natom)) !check deallocate
allocate(occupied(natom))
allocate(minimumIDs(natom))
allocate(unoccupied(natom))
allocate(CN_typeNo(natom,elemtype))
occupied = 0 ! keep global site ID

do 333 nat = 1,natom !over each atom

	unoccupiedCount = 0
	CN_typeNo = 0 ! number of different types
	!if(nat .eq.3 .or.nat .eq.4) print*,"*** atom ID: ",nat
	!if (Mod(nat,2000) .eq. 0) print*,"*** atom ID: ",nat

	do 221 i = 1,natom !loop over sites (not each atom)
	!if(nat .eq.3 .or.nat .eq.4)     print*," site i ",i,"occupied(i): ",occupied(i)
		if (occupied(i) .eq. 0)then ! site could insert the atom
			unoccupiedCount = unoccupiedCount + 1
	!	   print*," unoccupiedCount:  ",unoccupiedCount

			unoccupied(unoccupiedCount) = i
	!	   print*,"unocuupied(unoccupiedCount):  ",unoccupied(unoccupiedCount)

			do 222 j = 1,CN_No(i,1) !only consider the first neighbour atoms
	          JID = CN_ID(i,1,j) !the first neighbour sites of site i
	!if(nat .eq.3.or.nat .eq.4) 			print*," site j ",JID
	!if(nat .eq.3.or.nat .eq.4) 			print*," occupied(JID) ",occupied(JID)

				if (occupied(JID) .ne. 0)then ! The site has been occupied
					do k = 1,elemtype !classify the first neighbour atom types
	!if(nat .eq.3.or.nat .eq.4) 					print*," type k ",k

						if(temptype(JID).eq.k)then 
							!if(nat .eq.3.or.nat .eq.4) print*," unoccupiedCount,k: ",unoccupiedCount,k
							
							CN_typeNo(unoccupiedCount,k) = CN_typeNo(unoccupiedCount,k) + 1
							!if(nat .eq.3.or.nat .eq.4) print*," CN_typeNo(unoccupiedCount,k): ",CN_typeNo(unoccupiedCount,k)
						endif
					enddo
				endif
			222 continue

		endif	
	!if(nat .eq.3.or.nat .eq.4) print*,"loop done at sitenat: ", nat

	221  continue
	!if(nat .eq.3.or.nat .eq.4)print*,"***** atom ID: ",nat,"unoccupiedCount: ",unoccupiedCount," SUM: ",nat+unoccupiedCount-1

	itempType = atype(nat)
	itempid = minloc(CN_typeNo(1:unoccupiedCount,itempType),DIM = 1)! minimal the same type number
	!itempNo = minval(CN_typeNo(1:unoccupiedCount,itempType))! minimal the same type number
    !if(nat .eq.3.or.nat .eq.4)print *,"itempid",itempid,"unoccupied(itempid(1)): ",unoccupied(itempid(1))
	itempsite = unoccupied(itempid(1)) ! convert to real site ID
    !if(nat .eq.3.or.nat .eq.4)print *,"itempsite: ",itempsite
    !if(nat .eq.3.or.nat .eq.4)print *,"itempType: ",itempsite
	temptype(itempsite) = itempType
!	!if(nat .eq.3.or.nat .eq.4) print *, "Before occupied(itempsite) ",occupied(itempsite)
	occupied(itempsite) = 1 

	! find all IDs with minumal values
	!itempType = atype(nat)
	!itempNo = minval(CN_typeNo(1:unoccupiedCount,itempType))! minimal the same type number
	!minimumCounter = 0
	!do imini=1,unoccupiedCount
	!	if (CN_typeNo(imini,itempType) .eq. itempNo)then
	!		minimumCounter = minimumCounter + 1
	!		minimumIDs(minimumCounter) = imini ! keep all IDs with minimums
	!	endif
	!enddo
!!
	!call random_number(r)
  	!irand1 = int(r * minimumCounter) + 1 !select one with minimum (need be converted)
	!irand2 = minimumIDs(irand1)
	!itempsite = unoccupied(irand2) ! convert to real site ID
	!temptype(itempsite) = itempType
	!!print *, "Before occupied(itempsite) ",occupied(itempsite)
!!
	!occupied(itempsite) = 1 !occupied by an atom
	!print *, "itempsite: ",itempsite,"itempType: ",itempType
	!print *, "After occupied(itempsite): ",occupied(itempsite),"itempNo: ",itempNo
	!print *,""
	!print *,""
	!pause
333 continue ! insert atom type loop
atype = temptype

deallocate(temptype) 
deallocate(occupied)
deallocate(minimumIDs)
deallocate(unoccupied)
deallocate(CN_typeNo)

!call lmpdata("Rearrange",0)	

!print *, "Output rearrange data "
!pause

return 
end