#########################################################
# FreeBSD                                               #
#########################################################
OBJEXT=.o
LIBEXT=.a
EXEEXT= 
F2CEXT=.f
PATHSEP=/
DEFFLG=-D

FC        ?= f77
FFLAGS    += -Os -fno-second-underscore
FOUTFLG   =-o 

#CC        = cc
CFLAGS    += -Os -D_POSIX_C_SOURCE=199506L -fPIC
COUTFLG   = -o

LD        = $(CC) 
LDFLAGS   = $(CFLAGS) -static
LOUTFLG   = $(COUTFLG)

AR        = ar cr
#AOUTFLG   =

RANLIB    = ranlib
RM        = rm -rf

LIBBLAS   = -L/usr/local/lib -lf77blas -lcblas -latlas -lg2c
LIBLAPACK = -L/usr/local/lib -llapack

LIBMETIS  = -L/usr/local/lib -lmetis 

LIBF77 = -lg2c  
#compat is required for ftime()
LIBC   = -lm -lcompat 

#########################################################
