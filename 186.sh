#PBS -l nodes=1:ppn=1
#PBS -N equal_tension
##############################################################
#export PATH=/opt/mpich_download/mpich-3.3.1/mpich_install/bin:$PATH
#export LD_LIBRARY_PATH=/opt/mpich_download/mpich-3.3.1/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/opt/lammps_nnp_download/n2p2/lib:${LD_LIBRARY_PATH}
#
#cd $PBS_O_WORKDIR
#
#mpirun -np 32 /opt/lammps/lmp_mpi -in NPT2_modify.in
###########################################################################

export PATH=/opt/OpenMPI/OpenMPI-pgi/bin/:$PATH
export LD_LIBRARY_PATH=/opt/OpenMPI/OpenMPI-pgi/lib/:$LD_LIBRA 

cd $PBS_O_WORKDIR
perl 0l.pl
#/opt/lammps2/OpenMPI/lmp_pgi -sf omp -pk omp 16 -in NPT2_modify.in
#export PATH=/opt/mpich_download/mpich-3.3.1/mpich_install/bin:$PATH
#export LD_LIBRARY_PATH=/opt/mpich_download/mpich-3.3.1/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/opt/lammps_nnp_download/n2p2/lib:${LD_LIBRARY_PATH}


#mpirun -np 24 /opt/lammps/lmp_mpi -in Tension.in
#/opt/lammps2/OpenMPI/lmp_pgi -sf omp -pk omp 16 -in NPT2.in
#perl binary_npt.pl