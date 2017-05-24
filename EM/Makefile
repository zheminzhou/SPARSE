############
# Makefile #
############


#########################################################
# Please change values based on your system environment #
#########################################################

# C++ Compiler
CC = g++

# Open MPI Compiler
MPICC = mpiCC

# C++ Compile Flag
CFLAGS = -c -Wall

# GCC OpenMP Flag
OMPFLAG = -fopenmp

# Static build
STFLAG = -static

# Include directories 
INCL = -I./Ipopt/include

# Library
LDFLAGS = -L./Ipopt/lib64 -L/usr/lib64/gcc/x86_64-suse-linux/4.7 -L/usr/lib64/gcc/x86_64-suse-linux/4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib64/gcc/x86_64-suse-linux/4.7/../../../../x86_64-suse-linux/lib -L/usr/lib64/gcc/x86_64-suse-linux/4.7/../../.. -lipopt -llapack -ldl -lcoinmumps -lfblas -lpthread -lgfortran -lm -lquadmath


#########################################################
# Please do NOT change below Makefile variables         #
#########################################################

# Common source files
SIGMA_COMMON_SRC = utils.cpp \


# SIGMA solve model
SIGMA_SOLVE_MODEL_SRC = sigma_solve_model_main.cpp \
	sigma_ipopt_run.cpp \
	sigma_ipopt_nlp.cpp


# Object files are from source files
SIGMA_COMMON_OBJS = ${SIGMA_COMMON_SRC:.cpp=.o}
SIGMA_SOLVE_MODEL_OBJS = ${SIGMA_SOLVE_MODEL_SRC:.cpp=.o}
# This should be the name of your executable

SIGMA_SOLVE_MODEL_EXE = sigma-solve-model
# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

.SUFFIXES: .cpp .c .o .obj

.PHONY: all mpi

# Build arg

all: sigma_solve_model

sigma_solve_model: $(SIGMA_SOLVE_MODEL_EXE)

# Build

$(SIGMA_SOLVE_MODEL_EXE): $(SIGMA_COMMON_OBJS) $(SIGMA_SOLVE_MODEL_OBJS)
	$(CC) $(OMPFLAG) $(STFLAG) $^ $(LDFLAGS) -o $@


# Clean

clean:
	rm -f \
	$(SIGMA_SOLVE_MODEL_EXE) \
	$(SIGMA_COMMON_OBJS) \
	$(SIGMA_SOLVE_MODEL_OBJS) \
####
# To handle CC and MPICC, let compile individually
#.cpp.o:
#    $(CC) $(CFLAGS) $(EXTRA_FLAGS) $(INCL) -o $@ $<
#
#.cpp.obj:
#    $(CC) $(CFLAGS) $(EXTRA_FLAGS) $(INCL) $^ -o $@ `$(CYGPATH_W) '$<'`

utils.o: utils.cpp
	$(CC) $(CFLAGS) utils.cpp

sigma_solve_model_main.o: sigma_solve_model_main.cpp
	$(CC) $(CFLAGS) $(OMPFLAG) $(INCL) sigma_solve_model_main.cpp

sigma_ipopt_run.o: sigma_ipopt_run.cpp
	$(CC) $(CFLAGS) $(OMPFLAG) $(INCL) sigma_ipopt_run.cpp

sigma_ipopt_nlp.o: sigma_ipopt_nlp.cpp
	$(CC) $(CFLAGS) $(OMPFLAG) $(INCL) sigma_ipopt_nlp.cpp


