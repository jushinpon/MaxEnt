!aint(): truncates a real number by removing its fractional part
!The function NINT rounds to the nearest whole number    	
!The function INT will:Leave an integer unchanged.Round a real number towards 0.


!first, second, third, fourth neighbour IDs 
!CN_No(i,j)? i: atom id, j: neighbor type --> number of this neighbor type of atom i
!CN_ID(i,j,k)? i: atom id, j: neighbor type, k: counter --> ID
! allocate in information.f95
subroutine celllist
use information
implicit real*8(a-h,o-z)   
integer,allocatable::binpnt(:),bin(:)
integer nxcell,nycell,nzcell
dimension del(3)
! we need to modify the value of original rlist first
nxcell = 0
nycell = 0
nzcell = 0

do 21 i=1,natom
! all fractional coordinates should be positive
	 ix=dint(x(i)/rlist)+1
	 iy=dint(y(i)/rlist)+1
	 iz=dint(z(i)/rlist)+1
	 if(ix .ge. nxcell) nxcell = ix
	 if(iy .ge. nycell) nycell = iy
	 if(iz .ge. nzcell) nzcell = iz
21 continue

if(nxcell .lt. 3 .or. nycell .lt. 3 .or.nzcell .lt. 3) then
    write(*,*)"Cells are fewer than 3. STOP here!!"
	write(*,*)"nxcell",nxcell,"nycell",nycell,"nzcell",nzcell
	stop
endif

rlistsq = rlist*rlist

allocate(binpnt(nxcell*nycell*nzcell))
allocate(bin(natom))
       
binpnt=0 ! set the array to be 0 initially
bin =0

! need to do once only
do 2 i=1,natom
! all fractional coordinates should be positive
	 ix=aint(x(i)/rlist)+1
!if(i.eq.1 .or. i.eq. 401) write(*,*)"i= ",i," ix= ",ix,nxcell,x(i)
	 
	 if(ix .gt. nxcell)ix=nxcell
	 iy=aint(y(i)/rlist)+1
!	 if(i.eq.1 .or. i.eq. 401) write(*,*)"i= ",i," iy= ",iy,nycell
     if(iy .gt. nycell)iy=nycell 
	 iz=aint(z(i)/rlist)+1
!	 if(i.eq.1 .or. i.eq. 401) write(*,*)"i= ",i," iz= ",iz,nzcell
     if(iz .gt. nzcell)iz=nzcell

	 ib=(iz-1)*(nxcell*nycell)+(iy-1)*nxcell+ix
!.... ib is the cell ID atom i belongs to......
   	 bin(i)=binpnt(ib)
	 binpnt(ib)=i ! finally, binpnt(ib) keeps the lagest ID in this cell
!example, atoms 2,5,7 in cell 1 (ib = 1)
!ib=1,i = 2 --> bin(2) = binpnt(1) = 0 (first time), binpnt(1) = 2
!ib=1,i = 5 --> bin(5) = binpnt(1) = 2 (from the above value), binpnt(1) = 5
!ib=1,i = 7 --> bin(7) = binpnt(1) = 5 (from the above value), binpnt(1) = 7
! The final value of binpnt(1) after the natom loop will keep the larger ID in cell 1
2	continue

 CN_No = 0
 CN_ID = 0
	
	
! get the interaction energy here for Monte Carlo
do 3 i=1,natom
	 xtmp=x(i)													 
	 ytmp=y(i)
	 ztmp=z(i)
	
	 ixx=aint(xtmp/rlist)+1
	 if(ixx .gt. nxcell) ixx=nxcell !the atom is just located at the boundary
	 iyy=aint(ytmp/rlist)+1
	 if(iyy .gt. nycell)iyy=nycell
	 izz=aint(ztmp/rlist)+1
	 if(izz .gt. nzcell)izz=nzcell
!       write(*,*)"ixx,iyy,izz= ",ixx,iyy,izz
! the following goes through the j atom cell based on the reference atom i's cell
	  do 4 k1=-1,1
	    do 5 k2=-1,1
	      do 6 k3=-1,1
	       ix=ixx+k1
	       if(ix.lt.1) ix=nxcell
	       if(ix.gt.nxcell) ix=1

           iy=iyy+k2
	       if(iy.lt.1) iy=nycell
	       if(iy.gt.nycell) iy=1

	       iz=izz+k3
           if(iz.lt.1) iz=nzcell
	       if(iz.gt.nzcell) iz=1

     	   ib=(iz-1)*(nxcell)*(nycell)+(iy-1)*nxcell+ix !bin ID of atom j's cell
	       j=binpnt(ib)	! the higest ID in cell ib
!             write(*,*)"ib and j ",ib,j
10            if(j.ne.0)then  ! cell without atoms or atoms in this cell have been gone through        

                 if(j.ne.i)then !! i-j and j-i pairs are considered
                    dx=xtmp-x(j)
                    dy=ytmp-y(j)
                    dz=ztmp-z(j) 
           
      		   if (abs(dx) > half_xl .and. pbcx) then
                   dx = dx - sign(xl,dx)
               endif     
	           if (abs(dy) > half_yl .and. pbcy)then
                   dy = dy - sign(yl,dy)
               endif    
	           if (abs(dz) > half_zl .and. pbcz)then
                   dz = dz - sign(zl,dz)
               endif    
                    rsq=dx*dx + dy*dy + dz*dz
!write(*,*)"Check"
!write(*,*)i,j,rsq
                    if(rsq .le. rlistsq)then           

! begin classify the neighbor types of atom i                         index=index+1

                       r=dsqrt(rsq)
					   !write(*,*)i,j,r
                      ! do the classification for five types
					  vneimin = 1000.d0 ! initial value for finding the proper neighbor ID
					  neitypeID = 0 !integer
					  do icn =1,3
                        vnei = abs( ( r/rdfpeak(icn) ) - 1 )
						!write(*,*)"icn: ",icn
						!write(*,*)"vnei: ",vnei
						if(vnei .le. vneimin)then
						  vneimin = vnei
						  neitypeID = icn
						endif
					  enddo
					  if(neitypeID .eq. 0)then
					    write(*,*)"classifying neighbor type error!"
						write(*,*)"ref i: ",i
						write(*,*)"ref j: ",j
						write(*,*)"distance r: ",r
						!pause					    
					  else
					    CN_No(i,neitypeID) = CN_No(i,neitypeID) + 1
					    itemp = CN_No(i,neitypeID)
                        CN_ID(i,neitypeID,itemp) = j
					  endif
	                                         
                    endif !rsq .le. rlistsq
                  

                  endif !j.ne.i
                   j=bin(j)
				 goto 10

	            endif!j.ne.0
6	      continue ! loop of z cell
5	    continue ! loop of y cell
4	  continue ! loop of x cell

!      ncount(i)=index
3 continue ! loop of natom

deallocate(binpnt)
deallocate(bin)
!..................... finish cell list
! find the number of each second nearest atoms
!do i=1,10
!write(*,*)i,CN_No(i,1)
!write(*,*)i,CN_No(i,2)
!   do j = 1,CN_No(i,2)
!	write(*,*)"**neigh 2",j,": ",CN_ID(i,2,j)
!   enddo
!write(*,*)i,CN_No(i,3)
!write(*,*)i,CN_No(i,4)
!write(*,*)i,CN_No(i,5)
! 
!write(*,*)
!!if(i .eq. 2) then
!!	pause
!!endif
!enddo
!pause
return
end  
