CFLAGS	         =  
FFLAGS	         =-Wno-tabs
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = src/ksp/ksp/examples/tutorials/
MANSEC           = KSP
CLEANFILES       = main*.o *.mod 
NP               = 1
OBJ				 = dmc_type.o dmc.o 
OBJMAIN			 = ${OBJ} main.o 

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: ${OBJMAIN}  chkopts
	-${FLINKER} -o main ${OBJMAIN}  ${PETSC_KSP_LIB}

small:
	make clean
	make main 
	-@${MPIEXEC} -n 4 ./main

middle:
	make clean
	make main 
	-@${MPIEXEC} -n 4 ./main -log_view


#include ${PETSC_DIR}/lib/petsc/conf/test
