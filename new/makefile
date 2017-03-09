
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

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

OBJS_FUNC =dm_func_consts.o dm_func_ones.o

OBJS  = func_helper.o \
	petsc_helper.o function_helper.o \
	dispmodule.o dm_print.o dm_common.o \
	dm_data.o c_gettimeofday.o dm_tensor.o \
	dm_expr.o ${OBJS_FUNC} dm_basics.o ot_init.o

OBJ_MAIN=${OBJS} main.o

#--------------------------------------------------------------------------

all:
	@echo "\033[1;32m Cleaning... \033[0m"
	@cd build && make clean
	@echo "\033[1;32m Done. \033[0m"
	@echo "\033[1;32m Preprocessing files... \033[0m"
	@sh ./pre.sh
	@echo "\033[1;32m Done. \033[0m"
	@cd build && make main
	cp build/main ./

main: ${OBJ_MAIN} chkopts
	-${FLINKER} -o main ${OBJ_MAIN} ${PETSC_KSP_LIB} -lstdc++

#clean:
#	@cd build && make clean
#	mpifort -o main ${OBJ_MAIN} ${PETSC_KSP_LIB} -lstdc++

#include ${PETSC_DIR}/lib/petsc/conf/test

