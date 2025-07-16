#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=shared
#SBATCH --account=csd816
#SBATCH --tasks-per-node=8
#SBATCH -t 06:30:00
#SBATCH -J SHIRLEY_OCEAN



#OpenMP settings:
#OCEAN=/global/homes/s/sjamnuch/Perlmutter/OCEAN/bin/ocean.pl

#OpenMP settings:
ulimit -s unlimited
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=true
export SLURM_CPU_BIND="cores"

module purge
module load slurm/expanse/ DefaultModules sdsc/1.0 shared slurm/expanse/
module load cpu/0.15.4 gcc/10.2.0 openmpi/4.0.4 fftw/3.3.8
#module load cray-fftw
#module load espresso/7.0-libxc-5.2.2-cpu
#module list

#run the application:
#srun /global/u2/s/sjamnuch/Perlmutter/exciting/SHG/bin/excitingmpi > $SLURM_SUBMIT_DIR/OUT
#srun -n 64 $PARALLEL -in in.PTFE_C3F8x200 -screen none -log PTFE_C3F8x200.lammps.log
start=$(date +%s)
echo "Running scf"
srun -n 8 /home/sjamnuch/OCEAN/QE7/bin/pw.x  -inp scf.in >  scf.out
end=$(date +%s)
echo "Runtime: $((end - start )) sec"

echo "Running nscf"
start=$(date +%s)
srun -n 8 /home/sjamnuch/OCEAN/QE7/bin/pw.x  -inp nscf.in > nscf.out
end=$(date +%s)
echo "Runtime: $((end - start )) sec"

start=$(date +%s)
srun -n 8 /home/sjamnuch/OCEAN/QE7/bin/obf_basis.x < obf_basis.in > obf_basis.out
nbasis=`tail -n 6 obf_basis.out | head -n 1 | awk '{print $3}'`
echo "! ${sband} ${nbasis}"
end=$(date +%s)
echo "Runtime: $((end - start )) sec"

echo "Running ham"
start=$(date +%s)
srun -n 8 /home/sjamnuch/OCEAN/QE7/bin/obf_ham.x < obf_ham.in > obf_ham.out
end=$(date +%s)
echo "Runtime: $((end - start )) sec"


echo "Running diag at NBASIS: ${nbasis}" 
start=$(date +%s)
srun -n 8 /home/sjamnuch/OCEAN/QE7/bin/obf_diag.x  < obf_diag.in > obf_diag.out
end=$(date +%s)
echo "Runtime: $((end - start )) sec"

