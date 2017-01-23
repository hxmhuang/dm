# Distribted Matrix Wrapper 
It is a distributed matrix (DM) wrapper to simulate MATLAB environment on cluster platform.
This software is developed based on PETSc.

---

#Installation Guide:
## Pre-installation 
> Some third-party libraries must be installed including MPI compiler, PETSC and [pnetcdf](https://trac.mcs.anl.gov/projects/parallel-netcdf/wiki/Download) library, this installation only cover the process of PETSC and DM library installation.

## Install PETSC
1. Download PETSC using git
   
   `$ git clone -b maint https://bitbucket.org/petsc/petsc petsc`
	   
2. Set environmental variable before installing PETSC
   
   `$ export PETSC_DIR=[YOUR PETSC DIR]` </br>
   `$ export PETSC_ARCH=linux-gnu`
   
   Note: replace `[YOUR PETSC DIR]` with the full path of petsc directory created in previous step
3. Configure PETSC and make 
   
   `./configure --with-cc=mpicc --with-cxx=mpic++ --with-fc=mpifort --download-fblaslapack` </br>
   `make all test`

## Install DM library

1. Because DM library depends on pnetcdf, the path of pnetcdf library should be specified before compilation
   
	`$ export PATH_PNETCDF=[YOUR PNETCDF PATH]`
	
2. Download DM library from github
	
   `$ git clone https://github.com/hxmhuang/dm.git .`
   `$ git checkout 3D`	
   
3. There are many test cases in the DM source code, if you want to review the results of test cases, run 

	`$ cd dm`	
   `$ make small` 
   
   After the source code has been successfully compiled, you may see the output the test cases.
4. In order to use DM library in you own code, you must create a static DM library using command
   
   `$ make lib`
   
   If the compilation is successful, you can see the module files `*.mod` and a static library file
   `"libdm.a"` generated in DM directory. 

## Compile your own code that depend on DM library

Your have to specify the DM library path and PNETCDF path by setting the environment variable
 
>   `export PATH_DM=[YOUR DM LIBRARY]` </br>
>   `export PATH_PNETCDF=[YOUR PNETCDF PATH]`
   
   Here is an example of [makefile](./examples/makefile) you can use in your own project.

