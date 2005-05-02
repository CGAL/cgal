#########################################################
# Linux                                                 #
# Intel Compilers                                       #
# The C compiler defines __INTEL_COMPILER               #
#########################################################
OBJEXT=.o
LIBEXT=.a
EXEEXT= 
F2CEXT=.f
PATHSEP=/
DEFFLG=-D

FC        = ifc
FFLAGS    = -O3 
FOUTFLG   =-o ./

# -Xc: strict ANSI (Xa is extended, -c99 is C99)
# -c99: some c99 features (-c99- to disallow)
# -axW: generate Pentium4 optimized code, as well as generic 386
# -axK: generate Pentium3 optimized code, as well as generic 386
# -xK, -XW: generate Pentium? optimized code only
# -vec_report0: do not report information on vectorization
# -ansi_alias: assume no strange aliases (int to float, etc)
# -fno-fnalias: no array aliasing within functions
CC        = icc
CFLAGS    = -O3 -D_POSIX_C_SOURCE=199506L -c99
CFLAGS    = -O3 -D_POSIX_C_SOURCE=199506L -Xc -axW -ansi_alias -fno-fnalias -w1 -Werror
CFLAGS    = -c99 -O3 -D_POSIX_C_SOURCE=199506L -Xc \
            -xK \
            -vec_report0 \
            -ansi_alias -fno-fnalias \
            -w1 -Werror 
COUTFLG   = -o ./

LD        = $(CC) 
LDFLAGS   = $(CFLAGS)
LOUTFLG   = $(COUTFLG)

AR        = ar cr
AOUTFLG   =

RANLIB    = ranlib
RM        = rm -rf

LIBBLAS   = external/lib/linux/blas_aux.o \
            -L external/lib/linux -lf77blas -lcblas -latlas
LIBLAPACK = -L external/lib/linux -llapack
LIBMETIS  = -L external/lib/linux -lmetis 

LIBF77 = -lCEPCF90 -lIEPCF90 -lintrins -lF90 -limf -lpthread 
LIBC   = 

#########################################################







