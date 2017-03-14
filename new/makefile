
CFLAGS           =
FFLAGS		 = -fcheck=all -fbacktrace -J /Users/sxzhang/Projects/dm/new/build/
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

OBJS_FUNC = $(addprefix ./func/, consts.o ones.o rand.o ot_all_func.o)

OBJS_DATA = $(addprefix ./data/, data_consts.o data_rand.o )

OBJS_UTIL = $(addprefix ./utils/,func_ptr.o gettimeofday.o petsc_helper.o \
	    ot_display.o )

ABC = $(addprefix ./, ot_data.o  ot_tensor.o \
	ot_expr.o  ot_init.o)

OBJS  = ot_print.o ot_common.o ot_type.o ot_vector.o \
	${OBJS_UTIL} ${OBJS_DATA} ${ABC} ${OBJS_FUNC}$

OBJ_MAIN= ${OBJS} main.o

#--------------------------------------------------------------------------

all:
	@rm -rf main
	@echo ${OBJ_MAIN}
	@echo "\033[1;32m Cleaning... \033[0m"
	@cd build && make clean && rm -rf main
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

