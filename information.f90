! define all global parameters
!1.nx,ny,nz: replicated unit cell numbers in x, y, and z dimensions
!2. natom_unit: atom number in the unit cell (read input file)
!3.

MODULE INFORMATION
implicit real*8(a-h,o-z)
integer nx,ny,nz,natom_unit,nkMCshift ! replicate unit cell times in x, y, z dimensions,atom number in unit cell 
integer elemtype!total element type number, atom type
logical pbcx, pbcy, pbcz
real*8,allocatable::x(:),y(:),z(:),xkeep(:),ykeep(:),zkeep(:),atkeep(:),minatkeep(:)
real*8,allocatable::xmin(:),ymin(:),zmin(:),atomentropy(:),katomentropy(:),pairweight(:,:)
integer,allocatable::CN_No(:,:),CN_ID(:,:,:),atype(:)
integer natom,iterMC,inirand
character*30 name
real*8 rdfpeak(5),weight(5),ux(10),uy(10),uz(10) ! the first five rdf peaks from MS and entropy weight
real*8 xl,yl,zl,half_xl,half_yl,half_zl,rlist,confentropy,real_lattice
real*8 Entmin,kT,Entkeep,acceptratio,filterValue ! for MC
integer nadjust,naccept,keepNo ! for MC and counter for keeping ID with the large scoring values 
integer,allocatable:: keepeID(:) ! keep ID with the large scoring values

contains
!*********************************************
subroutine general_para
implicit real*8(a-h,o-z)
real shift
logical assigntype
real,allocatable::frac(:) ! atom fraction coordinate shift
character*128 buildWay,filename,tempChar 

call system('rm -rf output')
call system('mkdir output')
     
open(112,file="00input.dat",status='old')
read(112,*) !skip comment (# The input file for...)
read(112,*) !skip comment (####### begin below)
read(112,*) ! the total MC iteration times,the total initial random search times
read(112,*) iterMC,inirand
read(112,*) ! #! PBC in x, y, and z
read(112,*) pbcx,pbcy,pbcz
read(112,*)! rlist for cell list 
read(112,*)rlist
write(*,*)'rlist: ', rlist
read(112,*)!element type for HEA
read(112,*)elemtype
write(*,*)'element type for HEA: ',elemtype
read(112,*)!element fractions

allocate (frac(elemtype))
allocate (pairweight(elemtype,elemtype)) !! weighting for different pairs

total = 0.0
do 7 i=1,elemtype
	read(112,*)frac(i) ! fraction of an element type
	write(*,*)"fraction ",i,frac(i)
    total = total + frac(i)
7 continue
frac = frac/total !normalized

read(112,*)!need to assign atom types or not
read(112,*)assigntype
write(*,*)'assign atom types: ',assigntype

read(112,*)!#! how you build your structure (cfg, data, or default)
read(112,*)buildWay
write(*,*)'the way to build sructure: ',trim(buildWay)

read(112,*)!#! data file name, if buildWay is cfg or data
read(112,*)filename
write(*,*)'file name to read for cfg or data: ',trim(filename)

select case (trim(buildWay))
   
case ('cfg') 
      print*, "read cfg file"
    open(1111,file= trim(filename),status='old')  
	read(1111,*) !ITEM: TIMESTEP	
	read(1111,*) ! timestep
	read(1111,*) !ITEM: NUMBER OF ATOMS
	read(1111,*)natom
	read(1111,*) !ITEM: BOX BOUNDS pp pp pp
	read(1111,*)xlo, xhi
	read(1111,*)ylo, yhi
	read(1111,*)zlo, zhi
	read(1111,*) !ITEM: ATOMS id type x y z
	
close(112)  ! close the handle of 00input.dat 

case ('data')
      print*, "read data file" 
    open(111,file= trim(filename),status='old')  
	read(111,*) ! comment
	read(111,*)natom
	read(111,*) elemtype
	read(111,*)xlo, xhi
	read(111,*)ylo, yhi
	read(111,*)zlo, zhi
	
	print*,"xlo, xhi: ",xlo, xhi
	print*,"ylo, yhi: ",ylo, yhi
	print*,"zlo, zhi: ",zlo, zhi
	
close(112)  ! close the handle of 00input.dat      
      
      
case ('default') 
    print*, "use default way"
	read(112,*) ! atom number per unit cell
	read(112,*) natom_unit
	write(*,*)'! atom number per unit cell: ',natom_unit 
	read(112,*) ! fractional coordinates of each atom in the unit cell
	
	do 1 i=1,natom_unit
	read(112,*) ux(i),uy(i),uz(i) ! read the fractional coordinates of each atom in the unit cell
	write(*,*)"atom ID in unit cell: ",i
	1 continue
	read(112,*)! nx ny nz (the replicate times in x, y, and z dimensions)
	read(112,*)nx,ny,nz
	write(*,*)nx,ny,nz
	natom = natom_unit*nx*ny*nz         
	close(112)		
end select

! build the system

allocate (x(natom))
allocate (y(natom))
allocate (z(natom))

allocate (xkeep(natom))! configuration kept by MC
allocate (ykeep(natom))
allocate (zkeep(natom))

allocate (xmin(natom))
allocate (ymin(natom))
allocate (zmin(natom))

allocate (atype(natom))
allocate (atkeep(natom))
allocate (minatkeep(natom))

allocate (atomentropy(natom))
allocate (katomentropy(natom)) ! keep atom potential for the current best
allocate (keepeID(natom))

!!! get coordinates of all atoms

select case (trim(buildWay))
   
case ('cfg') 
	do 1001 i=1,natom
		read(1111,*)id_atom,itempatype,tempx,tempy,tempz !If id not sorted by lammps       
		x(id_atom)=tempx
		y(id_atom)=tempy
		z(id_atom)=tempz
		atype(id_atom)=itempatype
		write(*,*)id_atom,atype(id_atom),x(id_atom),y(id_atom),z(id_atom) !If id not sorted by lammps       
	1001 continue
    close(1111) ! close data file handle
!move all atoms within the box
	x = x - minval(x)
	y = y - minval(y)
	z = z - minval(z)
	
	xl = xhi-xlo
	yl = yhi-ylo
	zl = zhi-zlo
	half_xl = xl/2.d0    
	half_yl = yl/2.d0    
	half_zl = zl/2.d0
	
case ('data')
	read(111,*) tempChar!skip ATOMs
	!print*,tempChar
	do 100 i=1,natom
		read(111,*)id_atom,itempatype,tempx,tempy,tempz !If id not sorted by lammps       
		x(id_atom)=tempx
		y(id_atom)=tempy
		z(id_atom)=tempz
		atype(id_atom)=itempatype
		if(i == natom)then
			write(*,*)"show the coordinates of last atom for check"
			write(*,*)id_atom,atype(id_atom),x(id_atom),y(id_atom),z(id_atom) !If id not sorted by lammps       
        endif 
	100 continue
    close(111) ! close data file handle
!move all atoms within the box
	x = x - minval(x)
	y = y - minval(y)
	z = z - minval(z)
	
	xl = xhi-xlo
	yl = yhi-ylo
	zl = zhi-zlo
	half_xl = xl/2.d0    
	half_yl = yl/2.d0    
	half_zl = zl/2.d0
	
case ('default')
	IDcounter = 0
	do 2 i = 1,nx
		do 3 j = 1,ny
			do 4 k = 1,nz
				do 5 l =1,natom_unit
					IDcounter = IDcounter + 1
					x(IDcounter)= ux(l) + (i-1) 
					y(IDcounter)= uy(l) + (j-1) 
					z(IDcounter)= uz(l) + (k-1)			
				5 continue
			4 continue
		3 continue
	2 continue
	
	! the half box lengths in x, y, and z dimensions (based on fractional coordinates) 
	half_xl = dble(nx)/2.d0    
	half_yl = dble(ny)/2.d0    
	half_zl = dble(nz)/2.d0
	xl = dble(nx)
	yl = dble(ny)
	zl = dble(nz)
end select 

! assign the initial values for data output by lmpdata.f95
xkeep = x
ykeep = y
zkeep = z

xmin = x
ymin = y
zmin = z

atkeep =atype ! keep the atom types


!assign atom types
if(assigntype)then
	ncounter = 0 
	do 11 i=1,elemtype-1
	ielement = nint( dble(natom)*frac(i) ) ! round to the closest integer
	do j = 1,ielement
		ncounter = ncounter + 1 ! atom counter for assigned atom type
		atype(ncounter) = i !assign atom type
	enddo
	11 continue
	! assign the last element type
	do 12 i = ncounter+1,natom
		atype(i) = elemtype ! the last element type (largest type ID)
		!write(*,*)i,atype(i)
	12 continue
endif
	!write(*,*)ncounter, natom
	!pause  
! read unit cell atom fractional coordinates and nx, ny, nz

	  
allocate(CN_No(natom,5)) ! consider five neighbor types
allocate(CN_ID(natom,5,50)) ! consider five neighbor types and corresponding IDs.
! assume the maximal number is 50. 
CN_No = 0
CN_ID = 0
!write(*,*)'1'
! B2 rdf peaks by MS for a B2 unit cell (lattice constant = 1)
! 0.87 1.0 1.41 1.65 1.73 
!!!!!!!!!!!!!!
rdfpeak(1)=0.72
rdfpeak(2)=1.02
rdfpeak(3)=1.24
rdfpeak(4)=1.425
rdfpeak(5)=1.6
!write(*,*)'2'
weight = 1.d0
weightbase = 1.d0 
weight(1) = weightbase*1.e3 !first Id is the neighbour ID, the second is atom type
weight(2) = weightbase*0.0
weight(3) = weightbase*(-1.0)
weight(4) = weightbase*(-1.e-2) !not used 
weight(5) = weightbase*(-1.e-1) !not used

 !6 is the second nearest Number of a reference atom
pairweight = dble(weight(1)) ! initial values for all pairs
!pairweight = 1.0 ! initial values for all pairs
pairweight(1,2) = -4.566
pairweight(1,3) = -3.966
pairweight(1,4) = -5.467
pairweight(1,5) = -4.783
pairweight(2,1) = -4.566
pairweight(2,3) = -4.899
pairweight(2,4) = -6.373
pairweight(2,5) = -5.510
pairweight(3,1) = -3.966
pairweight(3,2) = -4.899
pairweight(3,4) = -4.985
pairweight(3,5) = -4.979
pairweight(4,1) = -5.467
pairweight(4,2) = -6.373
pairweight(4,3) = -4.985
pairweight(4,5) = -6.738
pairweight(5,1) = -4.783
pairweight(5,2) = -5.510
pairweight(5,3) = -4.979
pairweight(5,4) = -6.738
pairweight = pairweight/dble(weight(1))
!pairweight = 1.0

!write(*,*)'3'
!!!!!!!!!!!!!!

!write(*,*)'4' 	
call lmpdata("INI",0)	
endsubroutine general_para 
	 
ENDMODULE
