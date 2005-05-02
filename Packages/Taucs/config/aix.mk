#########################################################
# IBM POWER                                             #
#########################################################
CC        = xlc_r
CFLAGS    = -O3 -q64
FC        = xlf_r -O3 -q64
FFLAGS    = -O3 -q64
RM        = rm
AR        = ar -X 64 vq
RANLIB    = ranlib
LD        = xlc_r
LDFLAGS   = -g -q64 -bmaxdata:160000000000 -bmaxstack:160000000000
INCDIRS   = 
EXPFILES  = 
BLASLIB   = -lesslsmp
LAPACKLIB = 
OTHERLIBS = -lm_r -lxlf90_r
DEFINES   = -DNOMETIS
