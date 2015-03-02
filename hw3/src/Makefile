#===============================================================================
# User Options
#===============================================================================

COMPILER    = gnu
MPI         = no
OPENMP      = yes
OPTIMIZE    = yes
DEBUG       = no
PROFILE     = no
PAPI        = no
BENCHMARK   = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = diffSolver

source = \
main.cpp \
diffSolver.cpp \
Sparse.cpp \
Mesh.cpp \
XSdata.cpp \
utils.cpp

obj = $(source:.cpp=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = g++
endif

# Standard Flags
CFLAGS := -std=c++11

# Profliing Flags
ifeq ($(PROFILE),yes)
  CFLAGS += -pg
endif

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj)
	$(CC) $(obj) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj)

edit:
	vim -p $(source) *.h

run:
	./$(program)
