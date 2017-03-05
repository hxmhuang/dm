
CFLAGS           =
FFLAGS		 = -fcheck=all -fbacktrace
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

#--------------------------------------------------------------------------

all:
	@make clean
	echo "preprocessing files..."
	@sh ./pre.sh
	@cd build && make all
	cp build/main ./

clean:
	cd build && make clean
	cp build/main ./
#--------------------------------------------------------------------------

# PYTHON = python

# TESTEXAMPLES_C_PYTHON       = ex1.PETSc runex1 runex1py ex1.rm
# TESTEXAMPLES_FORTRAN_PYTHON = ex1f.PETSc runex1f ex1f.rm

#--------------------------------------------------------------------------

#include ${PETSC_DIR}/lib/petsc/conf/test
