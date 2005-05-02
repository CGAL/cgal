#########################################################
# IBM POWER                                             #
#########################################################
CC        = xlc_r
CFLAGS    = -g -D_POSIX_C_SOURCE=199506L
F77       = xlf
FFLAGS    = -g
EXESUFFIX =
CTIME     = utils/ibm_ftime_timer.c
SSOURCES  = utils/ibm_rs6000_timer.s
SSOURCES  = 
RM        = rm
AR        = ar vq
RANLIB    = echo
LD        = $(CC)
INCDIRS   = 
EXPFILES  = -bI:raw/s_io.exp
EXPFILES  = 
BLASLIB   = -L /home/c1majo/lib -lesslp2n
BLASLIB   = -lblas -llapack
BLASLIB   = -lblas 
LAPACKLIB = /afs/watson.ibm.com/agora/share/power/mathlib/lapack/lapack.a
LAPACKLIB = $(HOME)/lib/lapack_rs6k.a
LAPACKLIB = $(HOME)/tools/LAPACK/lapack_RS6K.a
OTHERLIBS = -lpthreads -lm -lc_r -lbsd -lxlf90
DEFINES   = -DUSE_PTHREADS -DUSE_ESSL

