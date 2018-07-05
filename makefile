# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# This file contains some definitions used to compile arpack++ examples
# with the g++ compiler under linux.

PLAT = linux
CPP = g++-mp-7
# CPP = gcc

ARPACKPP_DIR = ${HOME}/Packages/arpackpp
ARPACKPP_INC = ${ARPACKPP_DIR}/include

# Defining ARPACK, LAPACK, UMFPACK, SUPERLU, BLAS and FORTRAN libraries.
# See the arpack++ manual or the README file for directions on how to 
# obtain arpack, umfpack and SuperLU packages. 
# UMFPACK_LIB and SUPERLU_LIB must be declared only if umfpack and superlu 
# are going to be used. Some BLAS and LAPACK fortran routines are 
# distributed along with arpack fortran code, but the user should verify 
# if optimized versions of these libraries are available before installing 
# arpack. The fortran libraries described below are those required to link
# fortran and c++ code using gnu g++ and f77 compiler under linux.
# Other libraries should be defined if the user intends to compile
# arpack++ on another environment.


BLAS_LIB = -L/opt/local/lib/lapack/ -lblas
LAPACK_LIB = -L/opt/local/lib/lapack/ -llapack
# SUPERLU_LIB = -L/opt/local/lib -lsuperlu
# ARPACK_LIB = -L/opt/local/lib -larpack

# BLAS_LIB = -lblas
# LAPACK_LIB = -llapack
ARPACK_LIB = -L/opt/local/lib -larpack

CPP_WARNINGS = -fpermissive -Wall -ansi
# CPP_WARNINGS =  -pedantic-errors
CPP_DEBUG = -g
CPP_OPTIM = -O2 -std=c++11

CPP_INC = -I${ARPACKPP_INC}

CPP_FLAGS = ${CPP_DEBUG} ${CPP_WARNINGS} ${CPP_OPTIM}

ALL_LIBS = ${BLAS_LIB} ${LAPACK_LIB} ${ARPACK_LIB}
# Putting all libraries together.

# Specify the path of the libraries in your environment.

# ALL_LIBS = -l/opt/local/lib/lapack/libblas.dylib -l/opt/local/lib/lapack/liiblapack.dylib -l/opt/local/lib/libsuperlu.lib

OBJS = main.o
TARGET = tJsquare

all:${TARGET}
$(TARGET):${OBJS}
	${CPP} ${CPP_INC} ${CPP_FLAGS} ${ALL_LIBS} $^ -o $@

main.o:main.cpp
	${CPP} ${CPP_INC} ${CPP_FLAGS} ${ALL_LIBS} -c $< -o $@

clean:
	rm -f ${TARGET} *.o log *.dat
	rm -rf *.dSYM
