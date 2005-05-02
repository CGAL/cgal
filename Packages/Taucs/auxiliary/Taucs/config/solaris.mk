#########################################################
# Sun Solaris UltraSparc                                #
#########################################################
OBJEXT=.o
LIBEXT=.a
EXEEXT= 
F2CEXT=.f
PATHSEP=/
DEFFLG=-D

# For SMP's use:
#CC        = cc -lm -xO4 -xtarget=ultra -xarch=v8plusa -dalign -fsimple=2 -xautopar -xreduction
#CC        = cc -xO4 -xtarget=ultra -xarch=v8plusa -dalign -fsimple=2 \
#            -xdepend -stackvar -I$(HOME)/include

# these definitions are for the Sun compilers and libraries
FC        = f77
FFLAGS    = -xO3
FOUTFLG   = -o ./

CC        = cc
CFLAGS    = -Xc -xO3
COUTFLG   = -o ./

LD        = $(CC) 
LDFLAGS   = $(CFLAGS)
LOUTFLG   = $(COUTFLG)

AR        = ar vq
AOUTFLG   =

RANLIB    = ranlib
RM        = /bin/rm

LIBBLAS   = -L /opt/SUNWspro/lib -lsunperf -lsunmath -lfsu
LIBLAPACK = 
LIBMETIS  = 

LIBF77    = -lF77 -lM77 
LIBC      = -lsunmath -lm -lcx -lc -lmalloc

#########################################################
# and these are for GCC & ATLAS                         #
# using GCC and ATLAS may not necessarily reduce        #
# performance                                           #
#########################################################

FC        = g77
FFLAGS    = -O3 -g -fno-second-underscore -Wall
FOUTFLG   =-o 

CC        = gcc
CFLAGS    = -O3 -g -D_POSIX_C_SOURCE=199506L -Wall -pedantic -ansi -D_GNU_SOURCE 
COUTFLG   = -o

LD        = $(CC) 
LDFLAGS   = $(CFLAGS)
LOUTFLG   = $(COUTFLG)

AR        = ar cr
AOUTFLG   =

RANLIB    = ranlib
RM        = rm -rf

LIBBLAS   = -Lexternal/lib/solaris -lf77blas -lcblas -latlas -lg2c
LIBLAPACK = -Lexternal/lib/solaris -llapack -lg2c

LIBMETIS  = -Lexternal/lib/solaris -lmetis 

LIBF77 = -lg2c  
LIBC   = -lm 

#########################################################







