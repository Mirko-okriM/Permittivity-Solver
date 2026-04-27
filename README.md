# About

v1.1.2 (23.04.2025) written by Mirko Siegert and Nils Kerkmann

F-PECS (**F**requency-dependent **P**ermittivity **E**lectrical **C**onductivity **S**imulation) is a finite-volume method digital rock physics code for simulating effective permittivity and electrical conductivity of binary rock images.

For more information about F-PECS, as well as the verification and validation of the results, please refer to our publication:  
In preparation. Will be provided upon publication.

Feel free to contact us. We will be happy to help:
Mirko Siegert: [mirko.siegert@web.de](mailto:mirko.siegert@web.de)  
Nils Kerkmann: [nils.kerkmann@ieg.fraunhofer.de](mailto:nils.kerkmann@ieg.fraunhofer.de)

Testers and feedback are very welcome!

# F-PECS - How to Use

**Disclaimer:** The following instructions and settings are the ones we use. On other systems, different settings or additional adjustments may be necessary.

## 0\. Dependencies

| **Dependency**                  | **Version**        |
| ------------------------------- | ------------------ |
| Linux system                    | openSUSE Leap 15.0 |
| Bash                            | 4.4.23(1)-release  |
| SLURM                           | 18.08.6-2          |
| OpenMP                          | 201511             |
| gcc, g++, gfortran (SUSE Linux) | 7.3.1 20180323     |
| PETSc                           | 3.23.5             |

Specific versions of the packages installed with PETSc:

| **Package** | **Version** |
| ----------- | ----------- |
| MPI         | 4           |
| MPICH       | 4.3.0       |
| cmake       | 3.10.2      |
| bison       | 3.0         |

## 1\. Setting up PETSc

1. Follow the download instructions: <https://petsc.org/main/install/download/>

2. Change to the new directory (cd) and configure PETSc as follows:  
```
./configure PETSC_ARCH=linux-gnu-complex-mpi --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-mpich --with-batch --with-scalar-type=complex --download-fblaslapack=1 --with-x11=0 --with-windows-graphics=0}
```
3. Follow the instructions in the console. For more information, please refer to the official PETSc installation guide:  <https://petsc.org/main/install/install/>

## 2\. Setup a Simulation

1. Create a case folder and paste the "src"-folder and all files including the ".raw" model file in there.
2. Open the "settings.F90" file and change entries for:
  - _casePath_ to path to case folder
  - Optional: Any other variable (see the relevant comment for an explanation)
3. Open the provided "run.sh" file and change entries for:
  - _ntasks-per-node_ to number of cores on node
  - _PETSC_DIR_ to main PETSc directory
  - _PETSC_ARCH_ to PETSc configuration
  - Path to the case folder

## 3\. Run Simulation

1. Submit the "run.sh" file as a batch job.
2. Results can be viewed in the ".out" file.
3. Any errors can be viewed in the ".err" file.
