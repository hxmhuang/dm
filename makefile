CFLAGS	         =  
FFLAGS	         =-Wno-tabs
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = src/ksp/ksp/examples/tutorials/
MANSEC           = KSP
CLEANFILES       = main*.o *.mod 
NP               = 1
OBJ				 = dm_type.o dm.o 
OBJMAIN			 = ${OBJ} main.o 

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: ${OBJMAIN}  chkopts
	-${FLINKER} -o main ${OBJMAIN}  ${PETSC_KSP_LIB}

small:
	make clean
	make main 
	-@${MPIEXEC} -n 4 ./main -m 3 -n 2 -ep 3.1

middle:
	make clean
	make main 
	-@${MPIEXEC} -n 4 ./main -m 3 -n 2 -ep 3.1 -debug -log_view


#include ${PETSC_DIR}/lib/petsc/conf/test
