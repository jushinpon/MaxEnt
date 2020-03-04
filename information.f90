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
real,allocatable::frac(:) ! atom fraction coordinate shift
 
!  B2 (000) (0.5,0.5,0.5)
shift = 0.25
!******** PBC conditions for your system *********
     
open(112,file="00input.dat",status='old')
read(112,*) !skip comment (unit cell structure, like B2)
read(112,*) ! the total MC iteration times,the total initial random search times
read(112,*) iterMC,inirand
read(112,*) ! atom number per unit cell
read(112,*) natom_unit
write(*,*)'! atom number per unit cell: ',natom_unit 
read(112,*) ! fractional coordinates of each atom in the unit cell
shift = 0.25 ! shift atom for cell list
do 1 i=1,natom_unit
   read(112,*) ux(i),uy(i),uz(i) ! read the fractional coordinates of each atom in the unit cell
   !ux(i) = ux(i) + shift
   !uy(i) = uy(i) + shift
   !uz(i) = uz(i) + shift
   write(*,*)"atom ID in unit cell: ",i
1 continue
read(112,*)! nx ny nz (the replicate times in x, y, and z dimensions)
read(112,*)nx,ny,nz
write(*,*)nx,ny,nz
read(112,*)! rlist for cell list (0.87 1.0 1.4 1.65 1.73 for B2. maximal value+0.2 should be used)
read(112,*)rlist
write(*,*)'rlist: ', rlist
read(112,*)! real lattice constant
read(112,*)real_lattice
write(*,*)'real lattice constant',real_lattice
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
close(112)
! build the system
natom = natom_unit*nx*ny*nz
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
IDcounter = 0
do 2 i = 1,nx
	do 3 j = 1,ny
		do 4 k = 1,nz
			do 5 l =1,natom_unit
				IDcounter = IDcounter + 1
				!write(*,*)"Assign atom",IDcounter
				x(IDcounter)= ux(l) + (i-1) 
				y(IDcounter)= uy(l) + (j-1) 
				z(IDcounter)= uz(l) + (k-1)			
			5 continue
		4 continue
	3 continue
2 continue	

! assign the initial values for data output by lmpdata.f95
xkeep = x
ykeep = y
zkeep = z

xmin = x
ymin = y
zmin = z

atkeep =atype ! keep the atom types


!assign atom types
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
	!write(*,*)ncounter, natom
	!pause  
! read unit cell atom fractional coordinates and nx, ny, nz
pbcx = .True.
pbcy = .True.
pbcz = .True.
	  
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
weight(1) = weightbase*1.e6 !first Id is the neighbour ID, the second is atom type
weight(2) = weightbase*(0.0)
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

!write(*,*)'3'
!!!!!!!!!!!!!!
! the half box lengths in x, y, and z dimensions (based on fractional coordinates) 
half_xl = dble(nx)/2.d0    
half_yl = dble(ny)/2.d0    
half_zl = dble(nz)/2.d0
xl = dble(nx)
yl = dble(ny)
zl = dble(nz)
    
!write(*,*)'4' 	
call lmpdata("INI",0)	
endsubroutine general_para 
	 
ENDMODULE
