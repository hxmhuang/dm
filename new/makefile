
CFLAGS           =
FFLAGS		 =
CPPFLAGS         =
FPPFLAGS         =
LOCDIR		 = tutorials/python/
MANSEC           =
EXAMPLESC	 = 
EXAMPLESF	 =
EXAMPLESCH	 =
EXAMPLESFH       =
EXAMPLESMATLAB   =
DIRS		 =
CLEANFILES       =

#PETSC_ARCH=linux-gnu-c-debug
#PETSC_DIR=/Users/sxzhang/Projects/petsc

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

OBJS_FUNC =dm_func_consts.o dm_func_ones.o
OBJS=dm_data.o dm_print.o c_gettimeofday.o dm_common.o dm_tensor.o dm_expr.o ${OBJS_FUNC} dm_basics.o
OBJ_MAIN=${OBJS} main.o

#--------------------------------------------------------------------------

all: ${OBJ_MAIN} chkopts
	-${FLINKER} -o main ${OBJ_MAIN} ${PETSC_KSP_LIB}
	-rm *.o *.mod

#--------------------------------------------------------------------------

PYTHON = python

TESTEXAMPLES_C_PYTHON       = ex1.PETSc runex1 runex1py ex1.rm
TESTEXAMPLES_FORTRAN_PYTHON = ex1f.PETSc runex1f ex1f.rm

#--------------------------------------------------------------------------

include ${PETSC_DIR}/lib/petsc/conf/test
