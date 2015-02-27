#===============================================================================
# User Options
#===============================================================================

COMPILER    = gnu
MPI         = no
OPENMP      = yes
OPTIMIZE    = yes
DEBUG       = yes
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


#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj)
	$(CC) $(obj) -o $@ -pg

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -Wall -pg

clean:
	rm -rf $(program) $(obj)

edit:
	vim -p $(source) *.h

run:
	./$(program)
