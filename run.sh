#!/bin/bash

#SBATCH --nodes=1				# DO NOT CHANGE (not ready for multiple nodes)
#SBATCH --ntasks-per-node=24	# Change to number of cores of used node
#SBATCH --cpus-per-task=1		# DO NOT CHANGE (>1 = multiple command executions)
#SBATCH --output=%x-%j.out   	# Output file named <jobname>-<jobid>.out
#SBATCH --error=%x-%j.err    	# Error file named <jobname>-<jobid>.err

export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE			# DO NOT CHANGE
export PETSC_DIR=/Path/to/PETSc/main/folder/			# Change to main path of PETSc
export PETSC_ARCH=Name-of-PETSc-configuration			# Change to PETSc configuration name (MPI MUST BE ENABLED!)

echo "This script is running on: $(hostname -f)"
date

# Path to case folder
cd /Path/to/case/folder/

${PETSC_DIR}/${PETSC_ARCH}/bin/mpif90 -cpp -fopenmp -O3 \
	-I${PETSC_DIR}/include \
	-I${PETSC_DIR}/${PETSC_ARCH}/include \
	permittivitySolver.f90 \
	-L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc -lm \
	-Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib \
	-o solverRun

# DEBUG: Check MPI linking
#echo "Checking MPI libraries:"
#ldd ./solverRun | grep -i mpi

# Start program with the full path to mpirun
${PETSC_DIR}/${PETSC_ARCH}/bin/mpirun -np $SLURM_NTASKS_PER_NODE ./solverRun
        