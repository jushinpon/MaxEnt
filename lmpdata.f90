!INTENT(IN), INTENT(OUT) or INTENT(INOUT).
subroutine lmpdata(prefix,dataID)
use information
implicit real*8(a-h,o-z)
integer, INTENT(IN) :: dataID
character*20 fname
CHARACTER(LEN=*), INTENT(IN) ::	prefix

write(fname,'(a,a,i5.5,a)')prefix,'_',dataID,'.data'
open(112,file=fname,status='unknown')

write(112,*)'LAMMPS data file via lmpdata.f95'
write(112,*)' '
write(112,'(i8,1x,a)')natom,'atoms'
write(112,'(i8,1x,a)')elemtype,'atom types'
write(112,*)' '
write(112,'(f10.3,1x,f10.3,1x,a,1x,a)')0.0,dble(nx)*real_lattice,'xlo','xhi'
write(112,'(f10.3,1x,f10.3,1x,a,1x,a)')0.0,dble(ny)*real_lattice,'ylo','yhi'
write(112,'(f10.3,1x,f10.3,1x,a,1x,a)')0.0,dble(nz)*real_lattice,'zlo','zhi'
write(112,*)' '
write(112,*)'Atoms'
write(112,*)' '

do 1 i=1,natom
	write(112,'(i8,1x,i8,1x,f10.3,1x,f10.3,1x,f10.3)')i,atype(i),x(i)*real_lattice, &
                     y(i)*real_lattice,z(i)*real_lattice
	!write(*,*)i,real_lattice
    !write(*,*)x(i),y(i),z(i)	
1 continue
close(112)
return	
end