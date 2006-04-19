# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
include $(CGAL_MAKEFILE)
#include /home/mpouget/libs/CGAL-3.1/make/makefile_i686_Linux-2.6.12-10-386_g++-4.0.2

PROFOPT= -g
CFLAGS=${PROFOPT} 
#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
LAPACK_DIR=/usr/local/CLAPACK
LAPACK_INCDIRS = -I$(LAPACK_DIR)/SRC -I$(LAPACK_DIR)

CXXFLAGS = \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   $(LAPACK_INCDIRS) \
           ${CFLAGS} 

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#
F2CDIR  = $(LAPACK_DIR)/F2CLIBS
LAPACK_LDLIBS  = $(LAPACK_DIR)/lapack_LINUX.a \
          $(LAPACK_DIR)/blas_LINUX.a $(LAPACK_DIR)/tmglib_LINUX.a \
          $(F2CDIR)/libF77.a $(F2CDIR)/libI77.a -lm -lc
#GSL_LIBS=-L${GSL_DIR}/lib -lgsl -lgslcblas -lm

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) \
          ${LAPACK_LDLIBS} 


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                blind blind_1pt

BLIND_OBJS=options.o PolyhedralSurf.o blind.o
blind: $(BLIND_OBJS)
	$(CGAL_CXX) $(LIBPATH) -g -o blind.exe $(BLIND_OBJS) \
	$(LDFLAGS) 

blind_1pt: blind_1pt.o
	$(CGAL_CXX) $(LIBPATH) -o blind_1pt.exe blind_1pt.o \
	$(LDFLAGS) 

rmexe:
	rm *.exe

clean: \
        blind.clean blind_1pt.clean options.clean PolyhedralSurf.clean \
	rmexe


depend: 
	makedepend *.[Ch]

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

# DO NOT DELETE

blind_1pt.o: /usr/include/stdio.h /usr/include/features.h
blind_1pt.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
blind_1pt.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
blind_1pt.o: /usr/include/bits/typesizes.h /usr/include/libio.h
blind_1pt.o: /usr/include/_G_config.h /usr/include/wchar.h
blind_1pt.o: /usr/include/bits/wchar.h /usr/include/gconv.h
blind_1pt.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
blind_1pt.o: /usr/include/stdlib.h /usr/include/sys/types.h
blind_1pt.o: /usr/include/time.h /usr/include/endian.h
blind_1pt.o: /usr/include/bits/endian.h /usr/include/sys/select.h
blind_1pt.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
blind_1pt.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
blind_1pt.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
blind_1pt.o: /usr/include/alloca.h ../../include/CGAL/Monge_via_jet_fitting.h
blind_1pt.o: ../../include/CGAL/jet_fitting_3_assertions.h
blind_1pt.o: /usr/include/math.h /usr/include/bits/huge_val.h
blind_1pt.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
blind_1pt.o: LinAlg_lapack.h
blind.o: /usr/include/stdio.h /usr/include/features.h
blind.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
blind.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
blind.o: /usr/include/bits/typesizes.h /usr/include/libio.h
blind.o: /usr/include/_G_config.h /usr/include/wchar.h
blind.o: /usr/include/bits/wchar.h /usr/include/gconv.h
blind.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
blind.o: /usr/include/stdlib.h /usr/include/sys/types.h /usr/include/time.h
blind.o: /usr/include/endian.h /usr/include/bits/endian.h
blind.o: /usr/include/sys/select.h /usr/include/bits/select.h
blind.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
blind.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
blind.o: /usr/include/bits/sched.h /usr/include/alloca.h
blind.o: ../../include/CGAL/Monge_via_jet_fitting.h
blind.o: ../../include/CGAL/jet_fitting_3_assertions.h /usr/include/math.h
blind.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
blind.o: /usr/include/bits/mathcalls.h GSL.h /usr/include/gsl/gsl_linalg.h
blind.o: /usr/include/gsl/gsl_mode.h /usr/include/gsl/gsl_permutation.h
blind.o: /usr/include/gsl/gsl_types.h /usr/include/gsl/gsl_errno.h
blind.o: /usr/include/errno.h /usr/include/bits/errno.h
blind.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
blind.o: /usr/include/asm-i386/errno.h /usr/include/gsl/gsl_check_range.h
blind.o: /usr/include/gsl/gsl_vector.h
blind.o: /usr/include/gsl/gsl_vector_complex_long_double.h
blind.o: /usr/include/gsl/gsl_complex.h
blind.o: /usr/include/gsl/gsl_vector_long_double.h
blind.o: /usr/include/gsl/gsl_block_long_double.h
blind.o: /usr/include/gsl/gsl_vector_complex.h
blind.o: /usr/include/gsl/gsl_block_complex_long_double.h
blind.o: /usr/include/gsl/gsl_vector_complex_double.h
blind.o: /usr/include/gsl/gsl_vector_double.h
blind.o: /usr/include/gsl/gsl_block_double.h
blind.o: /usr/include/gsl/gsl_block_complex_double.h
blind.o: /usr/include/gsl/gsl_vector_complex_float.h
blind.o: /usr/include/gsl/gsl_vector_float.h
blind.o: /usr/include/gsl/gsl_block_float.h
blind.o: /usr/include/gsl/gsl_block_complex_float.h
blind.o: /usr/include/gsl/gsl_vector_ulong.h
blind.o: /usr/include/gsl/gsl_block_ulong.h
blind.o: /usr/include/gsl/gsl_vector_long.h /usr/include/gsl/gsl_block_long.h
blind.o: /usr/include/gsl/gsl_vector_uint.h /usr/include/gsl/gsl_block_uint.h
blind.o: /usr/include/gsl/gsl_vector_int.h /usr/include/gsl/gsl_block_int.h
blind.o: /usr/include/gsl/gsl_vector_ushort.h
blind.o: /usr/include/gsl/gsl_block_ushort.h
blind.o: /usr/include/gsl/gsl_vector_short.h
blind.o: /usr/include/gsl/gsl_block_short.h
blind.o: /usr/include/gsl/gsl_vector_uchar.h
blind.o: /usr/include/gsl/gsl_block_uchar.h
blind.o: /usr/include/gsl/gsl_vector_char.h /usr/include/gsl/gsl_block_char.h
blind.o: /usr/include/gsl/gsl_matrix.h
blind.o: /usr/include/gsl/gsl_matrix_complex_long_double.h
blind.o: /usr/include/gsl/gsl_matrix_complex_double.h
blind.o: /usr/include/gsl/gsl_matrix_complex_float.h
blind.o: /usr/include/gsl/gsl_matrix_long_double.h
blind.o: /usr/include/gsl/gsl_matrix_double.h
blind.o: /usr/include/gsl/gsl_matrix_float.h
blind.o: /usr/include/gsl/gsl_matrix_ulong.h
blind.o: /usr/include/gsl/gsl_matrix_long.h
blind.o: /usr/include/gsl/gsl_matrix_uint.h /usr/include/gsl/gsl_matrix_int.h
blind.o: /usr/include/gsl/gsl_matrix_ushort.h
blind.o: /usr/include/gsl/gsl_matrix_short.h
blind.o: /usr/include/gsl/gsl_matrix_uchar.h
blind.o: /usr/include/gsl/gsl_matrix_char.h /usr/include/gsl/gsl_eigen.h
blind.o: PolyhedralSurf.h PolyhedralSurf_operations.h PolyhedralSurf_rings.h
blind.o: options.h
GSL.o: /usr/include/gsl/gsl_linalg.h /usr/include/gsl/gsl_mode.h
GSL.o: /usr/include/gsl/gsl_permutation.h /usr/include/stdlib.h
GSL.o: /usr/include/features.h /usr/include/sys/cdefs.h
GSL.o: /usr/include/gnu/stubs.h /usr/include/sys/types.h
GSL.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
GSL.o: /usr/include/bits/typesizes.h /usr/include/time.h
GSL.o: /usr/include/endian.h /usr/include/bits/endian.h
GSL.o: /usr/include/sys/select.h /usr/include/bits/select.h
GSL.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
GSL.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
GSL.o: /usr/include/bits/sched.h /usr/include/alloca.h
GSL.o: /usr/include/gsl/gsl_types.h /usr/include/gsl/gsl_errno.h
GSL.o: /usr/include/stdio.h /usr/include/libio.h /usr/include/_G_config.h
GSL.o: /usr/include/wchar.h /usr/include/bits/wchar.h /usr/include/gconv.h
GSL.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
GSL.o: /usr/include/errno.h /usr/include/bits/errno.h
GSL.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
GSL.o: /usr/include/asm-i386/errno.h /usr/include/gsl/gsl_check_range.h
GSL.o: /usr/include/gsl/gsl_vector.h
GSL.o: /usr/include/gsl/gsl_vector_complex_long_double.h
GSL.o: /usr/include/gsl/gsl_complex.h
GSL.o: /usr/include/gsl/gsl_vector_long_double.h
GSL.o: /usr/include/gsl/gsl_block_long_double.h
GSL.o: /usr/include/gsl/gsl_vector_complex.h
GSL.o: /usr/include/gsl/gsl_block_complex_long_double.h
GSL.o: /usr/include/gsl/gsl_vector_complex_double.h
GSL.o: /usr/include/gsl/gsl_vector_double.h
GSL.o: /usr/include/gsl/gsl_block_double.h
GSL.o: /usr/include/gsl/gsl_block_complex_double.h
GSL.o: /usr/include/gsl/gsl_vector_complex_float.h
GSL.o: /usr/include/gsl/gsl_vector_float.h /usr/include/gsl/gsl_block_float.h
GSL.o: /usr/include/gsl/gsl_block_complex_float.h
GSL.o: /usr/include/gsl/gsl_vector_ulong.h /usr/include/gsl/gsl_block_ulong.h
GSL.o: /usr/include/gsl/gsl_vector_long.h /usr/include/gsl/gsl_block_long.h
GSL.o: /usr/include/gsl/gsl_vector_uint.h /usr/include/gsl/gsl_block_uint.h
GSL.o: /usr/include/gsl/gsl_vector_int.h /usr/include/gsl/gsl_block_int.h
GSL.o: /usr/include/gsl/gsl_vector_ushort.h
GSL.o: /usr/include/gsl/gsl_block_ushort.h
GSL.o: /usr/include/gsl/gsl_vector_short.h /usr/include/gsl/gsl_block_short.h
GSL.o: /usr/include/gsl/gsl_vector_uchar.h /usr/include/gsl/gsl_block_uchar.h
GSL.o: /usr/include/gsl/gsl_vector_char.h /usr/include/gsl/gsl_block_char.h
GSL.o: /usr/include/gsl/gsl_matrix.h
GSL.o: /usr/include/gsl/gsl_matrix_complex_long_double.h
GSL.o: /usr/include/gsl/gsl_matrix_complex_double.h
GSL.o: /usr/include/gsl/gsl_matrix_complex_float.h
GSL.o: /usr/include/gsl/gsl_matrix_long_double.h
GSL.o: /usr/include/gsl/gsl_matrix_double.h
GSL.o: /usr/include/gsl/gsl_matrix_float.h
GSL.o: /usr/include/gsl/gsl_matrix_ulong.h /usr/include/gsl/gsl_matrix_long.h
GSL.o: /usr/include/gsl/gsl_matrix_uint.h /usr/include/gsl/gsl_matrix_int.h
GSL.o: /usr/include/gsl/gsl_matrix_ushort.h
GSL.o: /usr/include/gsl/gsl_matrix_short.h
GSL.o: /usr/include/gsl/gsl_matrix_uchar.h /usr/include/gsl/gsl_matrix_char.h
GSL.o: /usr/include/gsl/gsl_eigen.h
LinAlg_lapack.o: /usr/include/stdlib.h /usr/include/features.h
LinAlg_lapack.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
LinAlg_lapack.o: /usr/include/sys/types.h /usr/include/bits/types.h
LinAlg_lapack.o: /usr/include/bits/wordsize.h /usr/include/bits/typesizes.h
LinAlg_lapack.o: /usr/include/time.h /usr/include/endian.h
LinAlg_lapack.o: /usr/include/bits/endian.h /usr/include/sys/select.h
LinAlg_lapack.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
LinAlg_lapack.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
LinAlg_lapack.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
LinAlg_lapack.o: /usr/include/alloca.h
options.o: /usr/include/ctype.h /usr/include/features.h
options.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
options.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
options.o: /usr/include/bits/typesizes.h /usr/include/endian.h
options.o: /usr/include/bits/endian.h /usr/include/string.h options.h
PolyhedralSurf.o: PolyhedralSurf.h /usr/include/stdlib.h
PolyhedralSurf.o: /usr/include/features.h /usr/include/sys/cdefs.h
PolyhedralSurf.o: /usr/include/gnu/stubs.h /usr/include/sys/types.h
PolyhedralSurf.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
PolyhedralSurf.o: /usr/include/bits/typesizes.h /usr/include/time.h
PolyhedralSurf.o: /usr/include/endian.h /usr/include/bits/endian.h
PolyhedralSurf.o: /usr/include/sys/select.h /usr/include/bits/select.h
PolyhedralSurf.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
PolyhedralSurf.o: /usr/include/sys/sysmacros.h
PolyhedralSurf.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
PolyhedralSurf.o: /usr/include/alloca.h /usr/include/stdio.h
PolyhedralSurf.o: /usr/include/libio.h /usr/include/_G_config.h
PolyhedralSurf.o: /usr/include/wchar.h /usr/include/bits/wchar.h
PolyhedralSurf.o: /usr/include/gconv.h /usr/include/bits/stdio_lim.h
PolyhedralSurf.o: /usr/include/bits/sys_errlist.h
PolyhedralSurf.o: /usr/include/stdlib.h /usr/include/features.h
PolyhedralSurf.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
PolyhedralSurf.o: /usr/include/sys/types.h /usr/include/bits/types.h
PolyhedralSurf.o: /usr/include/bits/wordsize.h /usr/include/bits/typesizes.h
PolyhedralSurf.o: /usr/include/time.h /usr/include/endian.h
PolyhedralSurf.o: /usr/include/bits/endian.h /usr/include/sys/select.h
PolyhedralSurf.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
PolyhedralSurf.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
PolyhedralSurf.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
PolyhedralSurf.o: /usr/include/alloca.h /usr/include/stdio.h
PolyhedralSurf.o: /usr/include/libio.h /usr/include/_G_config.h
PolyhedralSurf.o: /usr/include/wchar.h /usr/include/bits/wchar.h
PolyhedralSurf.o: /usr/include/gconv.h /usr/include/bits/stdio_lim.h
PolyhedralSurf.o: /usr/include/bits/sys_errlist.h
