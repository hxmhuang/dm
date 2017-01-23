
CFLAGS	         =  
FFLAGS	         =-Wno-tabs -I ${PATH_PNETCDF}/include -J ${DM} -fbacktrace
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = src/ksp/ksp/examples/tutorials/
MANSEC           = KSP
CLEANFILES       = main *.o
NP               = 1
OBJ		 = 

OBJMAIN		 = ${OBJ} FCLS_single_pixel.o
#

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main : ${OBJMAIN}  chkopts
	-${FLINKER} -o main ${OBJMAIN} ${PETSC_KSP_LIB} ${PATH_PNETCDF}/lib/libpnetcdf.a -J ${DM} ${DM}/libdm.a


#include ${PETSC_DIR}/lib/petsc/conf/test
