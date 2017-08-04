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
CFLAGS = -c -Wall -std=c++0x

# GCC OpenMP Flag
OMPFLAG = -fopenmp

# Static build
STFLAG = -static

# Include directories 
INCL = -I./Ipopt/include

# Library
# LDFLAGS = -L./Ipopt/lib64 -L/usr/lib64/gcc/x86_64-suse-linux/4.7 -L/usr/lib64/gcc/x86_64-suse-linux/4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib64/gcc/x86_64-suse-linux/4.7/../../../../x86_64-suse-linux/lib -L/usr/lib64/gcc/x86_64-suse-linux/4.7/../../.. -lipopt -llapack -ldl -lcoinmumps -lfblas -lpthread -lgfortran -lm -lquadmath
LDFLAGS = -L./Ipopt/lib64 -L/usr/lib/gcc/x86_64-linux-gnu/5 -L/usr/lib/ -L/lib/ -L/lib64/ -L/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -lipopt -llapack -ldl -lcoinmumps -lfblas -lpthread -lgfortran -lm -lquadmath


#########################################################
# Please do NOT change below Makefile variables         #
#########################################################

# solve model
SOLVE_MODEL_SRC = solve_model_main.cpp \
	ipopt_run.cpp \
	ipopt_nlp.cpp


# Object files are from source files
SOLVE_MODEL_OBJS = ${SOLVE_MODEL_SRC:.cpp=.o}
# This should be the name of your executable

SOLVE_MODEL_EXE = solve-model
# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

.SUFFIXES: .cpp .c .o .obj

.PHONY: all mpi

# Build arg

all: solve_model

solve_model: $(SOLVE_MODEL_EXE)

# Build

$(SOLVE_MODEL_EXE): $(SOLVE_MODEL_OBJS)
	$(CC) $(OMPFLAG) $(STFLAG) $^ $(LDFLAGS) -o $@


# Clean

clean:
	rm -f \
	$(SOLVE_MODEL_EXE) \
	$(SOLVE_MODEL_OBJS) \
####
# To handle CC and MPICC, let compile individually
#.cpp.o:
#    $(CC) $(CFLAGS) $(EXTRA_FLAGS) $(INCL) -o $@ $<
#
#.cpp.obj:
#    $(CC) $(CFLAGS) $(EXTRA_FLAGS) $(INCL) $^ -o $@ `$(CYGPATH_W) '$<'`

solve_model_main.o: solve_model_main.cpp
	$(CC) $(CFLAGS) $(OMPFLAG) $(INCL) solve_model_main.cpp

ipopt_run.o: ipopt_run.cpp
	$(CC) $(CFLAGS) $(OMPFLAG) $(INCL) ipopt_run.cpp

ipopt_nlp.o: ipopt_nlp.cpp
	$(CC) $(CFLAGS) $(OMPFLAG) $(INCL) ipopt_nlp.cpp


