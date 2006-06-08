# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
include $(CGAL_MAKEFILE)

PROFOPT= -g
CFLAGS=${PROFOPT} 
#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
LAPACK_INCDIRS = -I$(LAPACK_DIR)/SRC -I$(LAPACK_DIR)
GSL_INCDIR = -I$(GSL_DIR)/usr/local/gsl

CXXFLAGS = \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   $(LAPACK_INCDIRS) \
	   $(GSL_INCDIR) \
           ${CFLAGS} 

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#
F2CDIR  = $(LAPACK_DIR)/F2CLIBS
LAPACK_LDLIBS  = $(LAPACK_DIR)/lapack_LINUX.a \
          $(LAPACK_DIR)/blas_LINUX.a $(LAPACK_DIR)/tmglib_LINUX.a \
          $(F2CDIR)/libF77.a $(F2CDIR)/libI77.a -lm -lc
GSL_LIBS=-L${GSL_DIR}/lib -lgsl -lgslcblas -lm

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          ${LAPACK_LDLIBS} \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) \
          ${GSL_LIBS}	

#---------------------------------------------------------------------#

#                    target entries
#---------------------------------------------------------------------#

all:            \
                blind blind_1pt blind_gsl blind_1pt_gsl

BLIND_OBJS=options.o PolyhedralSurf.o blind.o
blind: $(BLIND_OBJS)
	$(CGAL_CXX) $(LIBPATH) -g -o blind.exe $(BLIND_OBJS) \
	$(LDFLAGS) 

blind_1pt: blind_1pt.o
	$(CGAL_CXX) $(LIBPATH) -o blind_1pt.exe blind_1pt.o \
	$(LDFLAGS) 

BLIND_OBJS_GSL=options.o PolyhedralSurf.o blind_gsl.o
blind_gsl: $(BLIND_OBJS_GSL)
	$(CGAL_CXX) $(LIBPATH) -g -o blind_gsl.exe $(BLIND_OBJS_GSL) \
	$(LDFLAGS) 

blind_1pt_gsl: blind_1pt_gsl.o
	$(CGAL_CXX) $(LIBPATH) -o blind_1pt_gsl.exe blind_1pt_gsl.o \
	$(LDFLAGS) 



rmexe:
	rm *.exe

clean: \
        blind.clean blind_1pt.clean options.clean PolyhedralSurf.clean \
	blind_gsl.clean blind_1pt_gsl.clean \
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
blind_1pt_gsl.o: /usr/include/stdio.h /usr/include/features.h
blind_1pt_gsl.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
blind_1pt_gsl.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
blind_1pt_gsl.o: /usr/include/bits/typesizes.h /usr/include/libio.h
blind_1pt_gsl.o: /usr/include/_G_config.h /usr/include/wchar.h
blind_1pt_gsl.o: /usr/include/bits/wchar.h /usr/include/gconv.h
blind_1pt_gsl.o: /usr/include/bits/stdio_lim.h
blind_1pt_gsl.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
blind_1pt_gsl.o: /usr/include/sys/types.h /usr/include/time.h
blind_1pt_gsl.o: /usr/include/endian.h /usr/include/bits/endian.h
blind_1pt_gsl.o: /usr/include/sys/select.h /usr/include/bits/select.h
blind_1pt_gsl.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
blind_1pt_gsl.o: /usr/include/sys/sysmacros.h
blind_1pt_gsl.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
blind_1pt_gsl.o: /usr/include/alloca.h
blind_1pt_gsl.o: ../../include/CGAL/Monge_via_jet_fitting_gsl.h
blind_1pt_gsl.o: ../../include/CGAL/jet_fitting_3_assertions.h
blind_1pt_gsl.o: /usr/include/math.h /usr/include/bits/huge_val.h
blind_1pt_gsl.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
blind_1pt_gsl.o: GSL.h /usr/include/gsl/gsl_linalg.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_mode.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_permutation.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_types.h /usr/include/gsl/gsl_errno.h
blind_1pt_gsl.o: /usr/include/errno.h /usr/include/bits/errno.h
blind_1pt_gsl.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
blind_1pt_gsl.o: /usr/include/asm-i386/errno.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_check_range.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector_complex_long_double.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_complex.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector_long_double.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_block_long_double.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector_complex.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_block_complex_long_double.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector_complex_double.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector_double.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_block_double.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_block_complex_double.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector_complex_float.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector_float.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_block_float.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_block_complex_float.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector_ulong.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_block_ulong.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector_long.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_block_long.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector_uint.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_block_uint.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector_int.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_block_int.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector_ushort.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_block_ushort.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector_short.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_block_short.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector_uchar.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_block_uchar.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_vector_char.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_block_char.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_matrix.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_matrix_complex_long_double.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_matrix_complex_double.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_matrix_complex_float.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_matrix_long_double.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_matrix_double.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_matrix_float.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_matrix_ulong.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_matrix_long.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_matrix_uint.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_matrix_int.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_matrix_ushort.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_matrix_short.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_matrix_uchar.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_matrix_char.h
blind_1pt_gsl.o: /usr/include/gsl/gsl_eigen.h
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
blind.o: /usr/include/boost/property_map.hpp /usr/include/boost/config.hpp
blind.o: /usr/include/boost/config/user.hpp
blind.o: /usr/include/boost/config/select_compiler_config.hpp
blind.o: /usr/include/boost/config/compiler/gcc.hpp
blind.o: /usr/include/boost/config/select_stdlib_config.hpp
blind.o: /usr/include/boost/config/select_platform_config.hpp
blind.o: /usr/include/boost/config/posix_features.hpp /usr/include/unistd.h
blind.o: /usr/include/bits/posix_opt.h /usr/include/bits/confname.h
blind.o: /usr/include/getopt.h /usr/include/boost/config/suffix.hpp
blind.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
blind.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
blind.o: /usr/include/bits/posix2_lim.h
blind.o: /usr/include/boost/pending/cstddef.hpp
blind.o: /usr/include/boost/detail/iterator.hpp /usr/include/boost/config.hpp
blind.o: /usr/include/boost/concept_check.hpp /usr/include/boost/iterator.hpp
blind.o: /usr/include/boost/type_traits/conversion_traits.hpp
blind.o: /usr/include/boost/type_traits/is_convertible.hpp
blind.o: /usr/include/boost/type_traits/detail/yes_no_type.hpp
blind.o: /usr/include/boost/type_traits/config.hpp
blind.o: /usr/include/boost/type_traits/is_array.hpp
blind.o: /usr/include/boost/type_traits/detail/bool_trait_def.hpp
blind.o: /usr/include/boost/type_traits/detail/template_arity_spec.hpp
blind.o: /usr/include/boost/mpl/int.hpp /usr/include/boost/mpl/int_fwd.hpp
blind.o: /usr/include/boost/mpl/aux_/adl_barrier.hpp
blind.o: /usr/include/boost/mpl/aux_/config/adl.hpp
blind.o: /usr/include/boost/mpl/aux_/config/msvc.hpp
blind.o: /usr/include/boost/mpl/aux_/config/intel.hpp
blind.o: /usr/include/boost/mpl/aux_/config/gcc.hpp
blind.o: /usr/include/boost/mpl/aux_/config/workaround.hpp
blind.o: /usr/include/boost/detail/workaround.hpp
blind.o: /usr/include/boost/mpl/aux_/nttp_decl.hpp
blind.o: /usr/include/boost/mpl/aux_/config/nttp.hpp
blind.o: /usr/include/boost/preprocessor/cat.hpp
blind.o: /usr/include/boost/preprocessor/config/config.hpp
blind.o: /usr/include/boost/mpl/aux_/integral_wrapper.hpp
blind.o: /usr/include/boost/mpl/integral_c_tag.hpp
blind.o: /usr/include/boost/mpl/aux_/config/static_constant.hpp
blind.o: /usr/include/boost/mpl/aux_/static_cast.hpp
blind.o: /usr/include/boost/mpl/aux_/template_arity_fwd.hpp
blind.o: /usr/include/boost/mpl/aux_/preprocessor/params.hpp
blind.o: /usr/include/boost/mpl/aux_/config/preprocessor.hpp
blind.o: /usr/include/boost/preprocessor/comma_if.hpp
blind.o: /usr/include/boost/preprocessor/punctuation/comma_if.hpp
blind.o: /usr/include/boost/preprocessor/control/if.hpp
blind.o: /usr/include/boost/preprocessor/control/iif.hpp
blind.o: /usr/include/boost/preprocessor/logical/bool.hpp
blind.o: /usr/include/boost/preprocessor/facilities/empty.hpp
blind.o: /usr/include/boost/preprocessor/punctuation/comma.hpp
blind.o: /usr/include/boost/preprocessor/repeat.hpp
blind.o: /usr/include/boost/preprocessor/repetition/repeat.hpp
blind.o: /usr/include/boost/preprocessor/debug/error.hpp
blind.o: /usr/include/boost/preprocessor/detail/auto_rec.hpp
blind.o: /usr/include/boost/preprocessor/tuple/eat.hpp
blind.o: /usr/include/boost/preprocessor/inc.hpp
blind.o: /usr/include/boost/preprocessor/arithmetic/inc.hpp
blind.o: /usr/include/boost/mpl/aux_/config/lambda.hpp
blind.o: /usr/include/boost/mpl/aux_/config/ttp.hpp
blind.o: /usr/include/boost/mpl/aux_/config/ctps.hpp
blind.o: /usr/include/boost/mpl/aux_/config/overload_resolution.hpp
blind.o: /usr/include/boost/type_traits/integral_constant.hpp
blind.o: /usr/include/boost/mpl/bool.hpp /usr/include/boost/mpl/bool_fwd.hpp
blind.o: /usr/include/boost/mpl/integral_c.hpp
blind.o: /usr/include/boost/mpl/integral_c_fwd.hpp
blind.o: /usr/include/boost/mpl/aux_/lambda_support.hpp
blind.o: /usr/include/boost/type_traits/detail/bool_trait_undef.hpp
blind.o: /usr/include/boost/type_traits/add_reference.hpp
blind.o: /usr/include/boost/type_traits/is_reference.hpp
blind.o: /usr/include/boost/type_traits/detail/type_trait_def.hpp
blind.o: /usr/include/boost/type_traits/detail/type_trait_undef.hpp
blind.o: /usr/include/boost/type_traits/ice.hpp
blind.o: /usr/include/boost/type_traits/detail/ice_or.hpp
blind.o: /usr/include/boost/type_traits/detail/ice_and.hpp
blind.o: /usr/include/boost/type_traits/detail/ice_not.hpp
blind.o: /usr/include/boost/type_traits/detail/ice_eq.hpp
blind.o: /usr/include/boost/type_traits/is_arithmetic.hpp
blind.o: /usr/include/boost/type_traits/is_integral.hpp
blind.o: /usr/include/boost/type_traits/is_float.hpp
blind.o: /usr/include/boost/type_traits/is_abstract.hpp
blind.o: /usr/include/boost/static_assert.hpp
blind.o: /usr/include/boost/type_traits/is_class.hpp
blind.o: /usr/include/boost/type_traits/is_union.hpp
blind.o: /usr/include/boost/type_traits/remove_cv.hpp
blind.o: /usr/include/boost/type_traits/broken_compiler_spec.hpp
blind.o: /usr/include/boost/type_traits/detail/cv_traits_impl.hpp
blind.o: /usr/include/boost/type_traits/intrinsics.hpp
blind.o: /usr/include/boost/mpl/identity.hpp
blind.o: /usr/include/boost/mpl/aux_/na_spec.hpp
blind.o: /usr/include/boost/mpl/lambda_fwd.hpp
blind.o: /usr/include/boost/mpl/void_fwd.hpp
blind.o: /usr/include/boost/mpl/aux_/na.hpp
blind.o: /usr/include/boost/mpl/aux_/na_fwd.hpp
blind.o: /usr/include/boost/mpl/aux_/lambda_arity_param.hpp
blind.o: /usr/include/boost/mpl/aux_/arity.hpp
blind.o: /usr/include/boost/mpl/aux_/config/dtp.hpp
blind.o: /usr/include/boost/mpl/aux_/preprocessor/enum.hpp
blind.o: /usr/include/boost/mpl/aux_/preprocessor/def_params_tail.hpp
blind.o: /usr/include/boost/mpl/limits/arity.hpp
blind.o: /usr/include/boost/preprocessor/logical/and.hpp
blind.o: /usr/include/boost/preprocessor/logical/bitand.hpp
blind.o: /usr/include/boost/preprocessor/identity.hpp
blind.o: /usr/include/boost/preprocessor/facilities/identity.hpp
blind.o: /usr/include/boost/preprocessor/empty.hpp
blind.o: /usr/include/boost/preprocessor/arithmetic/add.hpp
blind.o: /usr/include/boost/preprocessor/arithmetic/dec.hpp
blind.o: /usr/include/boost/preprocessor/control/while.hpp
blind.o: /usr/include/boost/preprocessor/list/fold_left.hpp
blind.o: /usr/include/boost/preprocessor/list/detail/fold_left.hpp
blind.o: /usr/include/boost/preprocessor/control/expr_iif.hpp
blind.o: /usr/include/boost/preprocessor/list/adt.hpp
blind.o: /usr/include/boost/preprocessor/detail/is_binary.hpp
blind.o: /usr/include/boost/preprocessor/detail/check.hpp
blind.o: /usr/include/boost/preprocessor/logical/compl.hpp
blind.o: /usr/include/boost/preprocessor/list/fold_right.hpp
blind.o: /usr/include/boost/preprocessor/list/detail/fold_right.hpp
blind.o: /usr/include/boost/preprocessor/list/reverse.hpp
blind.o: /usr/include/boost/preprocessor/control/detail/while.hpp
blind.o: /usr/include/boost/preprocessor/tuple/elem.hpp
blind.o: /usr/include/boost/preprocessor/arithmetic/sub.hpp
blind.o: /usr/include/boost/mpl/aux_/config/eti.hpp
blind.o: /usr/include/boost/concept_archetype.hpp
blind.o: ../../include/CGAL/Monge_via_jet_fitting.h
blind.o: ../../include/CGAL/jet_fitting_3_assertions.h /usr/include/math.h
blind.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
blind.o: /usr/include/bits/mathcalls.h LinAlg_lapack.h PolyhedralSurf.h
blind.o: PolyhedralSurf_operations.h PolyhedralSurf_rings.h options.h
blind_gsl.o: /usr/include/stdio.h /usr/include/features.h
blind_gsl.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
blind_gsl.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
blind_gsl.o: /usr/include/bits/typesizes.h /usr/include/libio.h
blind_gsl.o: /usr/include/_G_config.h /usr/include/wchar.h
blind_gsl.o: /usr/include/bits/wchar.h /usr/include/gconv.h
blind_gsl.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
blind_gsl.o: /usr/include/stdlib.h /usr/include/sys/types.h
blind_gsl.o: /usr/include/time.h /usr/include/endian.h
blind_gsl.o: /usr/include/bits/endian.h /usr/include/sys/select.h
blind_gsl.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
blind_gsl.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
blind_gsl.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
blind_gsl.o: /usr/include/alloca.h /usr/include/boost/property_map.hpp
blind_gsl.o: /usr/include/boost/config.hpp /usr/include/boost/config/user.hpp
blind_gsl.o: /usr/include/boost/config/select_compiler_config.hpp
blind_gsl.o: /usr/include/boost/config/compiler/gcc.hpp
blind_gsl.o: /usr/include/boost/config/select_stdlib_config.hpp
blind_gsl.o: /usr/include/boost/config/select_platform_config.hpp
blind_gsl.o: /usr/include/boost/config/posix_features.hpp
blind_gsl.o: /usr/include/unistd.h /usr/include/bits/posix_opt.h
blind_gsl.o: /usr/include/bits/confname.h /usr/include/getopt.h
blind_gsl.o: /usr/include/boost/config/suffix.hpp /usr/include/limits.h
blind_gsl.o: /usr/include/bits/posix1_lim.h /usr/include/bits/local_lim.h
blind_gsl.o: /usr/include/linux/limits.h /usr/include/bits/posix2_lim.h
blind_gsl.o: /usr/include/boost/pending/cstddef.hpp
blind_gsl.o: /usr/include/boost/detail/iterator.hpp
blind_gsl.o: /usr/include/boost/config.hpp
blind_gsl.o: /usr/include/boost/concept_check.hpp
blind_gsl.o: /usr/include/boost/iterator.hpp
blind_gsl.o: /usr/include/boost/type_traits/conversion_traits.hpp
blind_gsl.o: /usr/include/boost/type_traits/is_convertible.hpp
blind_gsl.o: /usr/include/boost/type_traits/detail/yes_no_type.hpp
blind_gsl.o: /usr/include/boost/type_traits/config.hpp
blind_gsl.o: /usr/include/boost/type_traits/is_array.hpp
blind_gsl.o: /usr/include/boost/type_traits/detail/bool_trait_def.hpp
blind_gsl.o: /usr/include/boost/type_traits/detail/template_arity_spec.hpp
blind_gsl.o: /usr/include/boost/mpl/int.hpp
blind_gsl.o: /usr/include/boost/mpl/int_fwd.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/adl_barrier.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/config/adl.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/config/msvc.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/config/intel.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/config/gcc.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/config/workaround.hpp
blind_gsl.o: /usr/include/boost/detail/workaround.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/nttp_decl.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/config/nttp.hpp
blind_gsl.o: /usr/include/boost/preprocessor/cat.hpp
blind_gsl.o: /usr/include/boost/preprocessor/config/config.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/integral_wrapper.hpp
blind_gsl.o: /usr/include/boost/mpl/integral_c_tag.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/config/static_constant.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/static_cast.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/template_arity_fwd.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/preprocessor/params.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/config/preprocessor.hpp
blind_gsl.o: /usr/include/boost/preprocessor/comma_if.hpp
blind_gsl.o: /usr/include/boost/preprocessor/punctuation/comma_if.hpp
blind_gsl.o: /usr/include/boost/preprocessor/control/if.hpp
blind_gsl.o: /usr/include/boost/preprocessor/control/iif.hpp
blind_gsl.o: /usr/include/boost/preprocessor/logical/bool.hpp
blind_gsl.o: /usr/include/boost/preprocessor/facilities/empty.hpp
blind_gsl.o: /usr/include/boost/preprocessor/punctuation/comma.hpp
blind_gsl.o: /usr/include/boost/preprocessor/repeat.hpp
blind_gsl.o: /usr/include/boost/preprocessor/repetition/repeat.hpp
blind_gsl.o: /usr/include/boost/preprocessor/debug/error.hpp
blind_gsl.o: /usr/include/boost/preprocessor/detail/auto_rec.hpp
blind_gsl.o: /usr/include/boost/preprocessor/tuple/eat.hpp
blind_gsl.o: /usr/include/boost/preprocessor/inc.hpp
blind_gsl.o: /usr/include/boost/preprocessor/arithmetic/inc.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/config/lambda.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/config/ttp.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/config/ctps.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/config/overload_resolution.hpp
blind_gsl.o: /usr/include/boost/type_traits/integral_constant.hpp
blind_gsl.o: /usr/include/boost/mpl/bool.hpp
blind_gsl.o: /usr/include/boost/mpl/bool_fwd.hpp
blind_gsl.o: /usr/include/boost/mpl/integral_c.hpp
blind_gsl.o: /usr/include/boost/mpl/integral_c_fwd.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/lambda_support.hpp
blind_gsl.o: /usr/include/boost/type_traits/detail/bool_trait_undef.hpp
blind_gsl.o: /usr/include/boost/type_traits/add_reference.hpp
blind_gsl.o: /usr/include/boost/type_traits/is_reference.hpp
blind_gsl.o: /usr/include/boost/type_traits/detail/type_trait_def.hpp
blind_gsl.o: /usr/include/boost/type_traits/detail/type_trait_undef.hpp
blind_gsl.o: /usr/include/boost/type_traits/ice.hpp
blind_gsl.o: /usr/include/boost/type_traits/detail/ice_or.hpp
blind_gsl.o: /usr/include/boost/type_traits/detail/ice_and.hpp
blind_gsl.o: /usr/include/boost/type_traits/detail/ice_not.hpp
blind_gsl.o: /usr/include/boost/type_traits/detail/ice_eq.hpp
blind_gsl.o: /usr/include/boost/type_traits/is_arithmetic.hpp
blind_gsl.o: /usr/include/boost/type_traits/is_integral.hpp
blind_gsl.o: /usr/include/boost/type_traits/is_float.hpp
blind_gsl.o: /usr/include/boost/type_traits/is_abstract.hpp
blind_gsl.o: /usr/include/boost/static_assert.hpp
blind_gsl.o: /usr/include/boost/type_traits/is_class.hpp
blind_gsl.o: /usr/include/boost/type_traits/is_union.hpp
blind_gsl.o: /usr/include/boost/type_traits/remove_cv.hpp
blind_gsl.o: /usr/include/boost/type_traits/broken_compiler_spec.hpp
blind_gsl.o: /usr/include/boost/type_traits/detail/cv_traits_impl.hpp
blind_gsl.o: /usr/include/boost/type_traits/intrinsics.hpp
blind_gsl.o: /usr/include/boost/mpl/identity.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/na_spec.hpp
blind_gsl.o: /usr/include/boost/mpl/lambda_fwd.hpp
blind_gsl.o: /usr/include/boost/mpl/void_fwd.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/na.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/na_fwd.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/lambda_arity_param.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/arity.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/config/dtp.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/preprocessor/enum.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/preprocessor/def_params_tail.hpp
blind_gsl.o: /usr/include/boost/mpl/limits/arity.hpp
blind_gsl.o: /usr/include/boost/preprocessor/logical/and.hpp
blind_gsl.o: /usr/include/boost/preprocessor/logical/bitand.hpp
blind_gsl.o: /usr/include/boost/preprocessor/identity.hpp
blind_gsl.o: /usr/include/boost/preprocessor/facilities/identity.hpp
blind_gsl.o: /usr/include/boost/preprocessor/empty.hpp
blind_gsl.o: /usr/include/boost/preprocessor/arithmetic/add.hpp
blind_gsl.o: /usr/include/boost/preprocessor/arithmetic/dec.hpp
blind_gsl.o: /usr/include/boost/preprocessor/control/while.hpp
blind_gsl.o: /usr/include/boost/preprocessor/list/fold_left.hpp
blind_gsl.o: /usr/include/boost/preprocessor/list/detail/fold_left.hpp
blind_gsl.o: /usr/include/boost/preprocessor/control/expr_iif.hpp
blind_gsl.o: /usr/include/boost/preprocessor/list/adt.hpp
blind_gsl.o: /usr/include/boost/preprocessor/detail/is_binary.hpp
blind_gsl.o: /usr/include/boost/preprocessor/detail/check.hpp
blind_gsl.o: /usr/include/boost/preprocessor/logical/compl.hpp
blind_gsl.o: /usr/include/boost/preprocessor/list/fold_right.hpp
blind_gsl.o: /usr/include/boost/preprocessor/list/detail/fold_right.hpp
blind_gsl.o: /usr/include/boost/preprocessor/list/reverse.hpp
blind_gsl.o: /usr/include/boost/preprocessor/control/detail/while.hpp
blind_gsl.o: /usr/include/boost/preprocessor/tuple/elem.hpp
blind_gsl.o: /usr/include/boost/preprocessor/arithmetic/sub.hpp
blind_gsl.o: /usr/include/boost/mpl/aux_/config/eti.hpp
blind_gsl.o: /usr/include/boost/concept_archetype.hpp
blind_gsl.o: ../../include/CGAL/Monge_via_jet_fitting_gsl.h
blind_gsl.o: ../../include/CGAL/jet_fitting_3_assertions.h
blind_gsl.o: /usr/include/math.h /usr/include/bits/huge_val.h
blind_gsl.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h GSL.h
blind_gsl.o: /usr/include/gsl/gsl_linalg.h /usr/include/gsl/gsl_mode.h
blind_gsl.o: /usr/include/gsl/gsl_permutation.h /usr/include/gsl/gsl_types.h
blind_gsl.o: /usr/include/gsl/gsl_errno.h /usr/include/errno.h
blind_gsl.o: /usr/include/bits/errno.h /usr/include/linux/errno.h
blind_gsl.o: /usr/include/asm/errno.h /usr/include/asm-i386/errno.h
blind_gsl.o: /usr/include/gsl/gsl_check_range.h /usr/include/gsl/gsl_vector.h
blind_gsl.o: /usr/include/gsl/gsl_vector_complex_long_double.h
blind_gsl.o: /usr/include/gsl/gsl_complex.h
blind_gsl.o: /usr/include/gsl/gsl_vector_long_double.h
blind_gsl.o: /usr/include/gsl/gsl_block_long_double.h
blind_gsl.o: /usr/include/gsl/gsl_vector_complex.h
blind_gsl.o: /usr/include/gsl/gsl_block_complex_long_double.h
blind_gsl.o: /usr/include/gsl/gsl_vector_complex_double.h
blind_gsl.o: /usr/include/gsl/gsl_vector_double.h
blind_gsl.o: /usr/include/gsl/gsl_block_double.h
blind_gsl.o: /usr/include/gsl/gsl_block_complex_double.h
blind_gsl.o: /usr/include/gsl/gsl_vector_complex_float.h
blind_gsl.o: /usr/include/gsl/gsl_vector_float.h
blind_gsl.o: /usr/include/gsl/gsl_block_float.h
blind_gsl.o: /usr/include/gsl/gsl_block_complex_float.h
blind_gsl.o: /usr/include/gsl/gsl_vector_ulong.h
blind_gsl.o: /usr/include/gsl/gsl_block_ulong.h
blind_gsl.o: /usr/include/gsl/gsl_vector_long.h
blind_gsl.o: /usr/include/gsl/gsl_block_long.h
blind_gsl.o: /usr/include/gsl/gsl_vector_uint.h
blind_gsl.o: /usr/include/gsl/gsl_block_uint.h
blind_gsl.o: /usr/include/gsl/gsl_vector_int.h
blind_gsl.o: /usr/include/gsl/gsl_block_int.h
blind_gsl.o: /usr/include/gsl/gsl_vector_ushort.h
blind_gsl.o: /usr/include/gsl/gsl_block_ushort.h
blind_gsl.o: /usr/include/gsl/gsl_vector_short.h
blind_gsl.o: /usr/include/gsl/gsl_block_short.h
blind_gsl.o: /usr/include/gsl/gsl_vector_uchar.h
blind_gsl.o: /usr/include/gsl/gsl_block_uchar.h
blind_gsl.o: /usr/include/gsl/gsl_vector_char.h
blind_gsl.o: /usr/include/gsl/gsl_block_char.h /usr/include/gsl/gsl_matrix.h
blind_gsl.o: /usr/include/gsl/gsl_matrix_complex_long_double.h
blind_gsl.o: /usr/include/gsl/gsl_matrix_complex_double.h
blind_gsl.o: /usr/include/gsl/gsl_matrix_complex_float.h
blind_gsl.o: /usr/include/gsl/gsl_matrix_long_double.h
blind_gsl.o: /usr/include/gsl/gsl_matrix_double.h
blind_gsl.o: /usr/include/gsl/gsl_matrix_float.h
blind_gsl.o: /usr/include/gsl/gsl_matrix_ulong.h
blind_gsl.o: /usr/include/gsl/gsl_matrix_long.h
blind_gsl.o: /usr/include/gsl/gsl_matrix_uint.h
blind_gsl.o: /usr/include/gsl/gsl_matrix_int.h
blind_gsl.o: /usr/include/gsl/gsl_matrix_ushort.h
blind_gsl.o: /usr/include/gsl/gsl_matrix_short.h
blind_gsl.o: /usr/include/gsl/gsl_matrix_uchar.h
blind_gsl.o: /usr/include/gsl/gsl_matrix_char.h /usr/include/gsl/gsl_eigen.h
blind_gsl.o: PolyhedralSurf.h PolyhedralSurf_operations.h
blind_gsl.o: PolyhedralSurf_rings.h options.h
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
PolyhedralSurf.o: PolyhedralSurf.h /usr/include/boost/property_map.hpp
PolyhedralSurf.o: /usr/include/boost/config.hpp
PolyhedralSurf.o: /usr/include/boost/config/user.hpp
PolyhedralSurf.o: /usr/include/boost/config/select_compiler_config.hpp
PolyhedralSurf.o: /usr/include/boost/config/compiler/gcc.hpp
PolyhedralSurf.o: /usr/include/boost/config/select_stdlib_config.hpp
PolyhedralSurf.o: /usr/include/boost/config/select_platform_config.hpp
PolyhedralSurf.o: /usr/include/boost/config/posix_features.hpp
PolyhedralSurf.o: /usr/include/unistd.h /usr/include/features.h
PolyhedralSurf.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
PolyhedralSurf.o: /usr/include/bits/posix_opt.h /usr/include/bits/types.h
PolyhedralSurf.o: /usr/include/bits/wordsize.h /usr/include/bits/typesizes.h
PolyhedralSurf.o: /usr/include/bits/confname.h /usr/include/getopt.h
PolyhedralSurf.o: /usr/include/boost/config/suffix.hpp /usr/include/limits.h
PolyhedralSurf.o: /usr/include/bits/posix1_lim.h
PolyhedralSurf.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
PolyhedralSurf.o: /usr/include/bits/posix2_lim.h
PolyhedralSurf.o: /usr/include/boost/pending/cstddef.hpp
PolyhedralSurf.o: /usr/include/boost/detail/iterator.hpp
PolyhedralSurf.o: /usr/include/boost/config.hpp
PolyhedralSurf.o: /usr/include/boost/concept_check.hpp
PolyhedralSurf.o: /usr/include/boost/iterator.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/conversion_traits.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_convertible.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/yes_no_type.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/config.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_array.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/bool_trait_def.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/template_arity_spec.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/int.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/int_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/adl_barrier.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/adl.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/msvc.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/intel.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/gcc.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/workaround.hpp
PolyhedralSurf.o: /usr/include/boost/detail/workaround.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/nttp_decl.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/nttp.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/cat.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/config/config.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/integral_wrapper.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/integral_c_tag.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/static_constant.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/static_cast.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/template_arity_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/preprocessor/params.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/preprocessor.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/comma_if.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/punctuation/comma_if.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/control/if.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/control/iif.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/logical/bool.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/facilities/empty.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/punctuation/comma.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/repeat.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/repetition/repeat.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/debug/error.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/detail/auto_rec.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/tuple/eat.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/inc.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/arithmetic/inc.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/lambda.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/ttp.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/ctps.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/overload_resolution.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/integral_constant.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/bool.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/bool_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/integral_c.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/integral_c_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/lambda_support.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/bool_trait_undef.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/add_reference.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_reference.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/type_trait_def.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/type_trait_undef.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/ice.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/ice_or.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/ice_and.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/ice_not.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/ice_eq.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_arithmetic.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_integral.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_float.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_abstract.hpp
PolyhedralSurf.o: /usr/include/boost/static_assert.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_class.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_union.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/remove_cv.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/broken_compiler_spec.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/cv_traits_impl.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/intrinsics.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/identity.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/na_spec.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/lambda_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/void_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/na.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/na_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/lambda_arity_param.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/arity.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/dtp.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/preprocessor/enum.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/preprocessor/def_params_tail.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/limits/arity.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/logical/and.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/logical/bitand.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/identity.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/facilities/identity.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/empty.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/arithmetic/add.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/arithmetic/dec.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/control/while.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/list/fold_left.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/list/detail/fold_left.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/control/expr_iif.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/list/adt.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/detail/is_binary.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/detail/check.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/logical/compl.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/list/fold_right.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/list/detail/fold_right.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/list/reverse.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/control/detail/while.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/tuple/elem.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/arithmetic/sub.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/eti.hpp
PolyhedralSurf.o: /usr/include/boost/concept_archetype.hpp
PolyhedralSurf.o: /usr/include/stdlib.h /usr/include/sys/types.h
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
PolyhedralSurf.o: /usr/include/boost/property_map.hpp
PolyhedralSurf.o: /usr/include/boost/config.hpp
PolyhedralSurf.o: /usr/include/boost/config/user.hpp
PolyhedralSurf.o: /usr/include/boost/config/select_compiler_config.hpp
PolyhedralSurf.o: /usr/include/boost/config/compiler/gcc.hpp
PolyhedralSurf.o: /usr/include/boost/config/select_stdlib_config.hpp
PolyhedralSurf.o: /usr/include/boost/config/select_platform_config.hpp
PolyhedralSurf.o: /usr/include/boost/config/posix_features.hpp
PolyhedralSurf.o: /usr/include/unistd.h /usr/include/features.h
PolyhedralSurf.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
PolyhedralSurf.o: /usr/include/bits/posix_opt.h /usr/include/bits/types.h
PolyhedralSurf.o: /usr/include/bits/wordsize.h /usr/include/bits/typesizes.h
PolyhedralSurf.o: /usr/include/bits/confname.h /usr/include/getopt.h
PolyhedralSurf.o: /usr/include/boost/config/suffix.hpp /usr/include/limits.h
PolyhedralSurf.o: /usr/include/bits/posix1_lim.h
PolyhedralSurf.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
PolyhedralSurf.o: /usr/include/bits/posix2_lim.h
PolyhedralSurf.o: /usr/include/boost/pending/cstddef.hpp
PolyhedralSurf.o: /usr/include/boost/detail/iterator.hpp
PolyhedralSurf.o: /usr/include/boost/config.hpp
PolyhedralSurf.o: /usr/include/boost/concept_check.hpp
PolyhedralSurf.o: /usr/include/boost/iterator.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/conversion_traits.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_convertible.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/yes_no_type.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/config.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_array.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/bool_trait_def.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/template_arity_spec.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/int.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/int_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/adl_barrier.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/adl.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/msvc.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/intel.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/gcc.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/workaround.hpp
PolyhedralSurf.o: /usr/include/boost/detail/workaround.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/nttp_decl.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/nttp.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/cat.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/config/config.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/integral_wrapper.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/integral_c_tag.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/static_constant.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/static_cast.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/template_arity_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/preprocessor/params.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/preprocessor.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/comma_if.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/punctuation/comma_if.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/control/if.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/control/iif.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/logical/bool.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/facilities/empty.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/punctuation/comma.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/repeat.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/repetition/repeat.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/debug/error.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/detail/auto_rec.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/tuple/eat.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/inc.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/arithmetic/inc.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/lambda.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/ttp.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/ctps.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/overload_resolution.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/integral_constant.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/bool.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/bool_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/integral_c.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/integral_c_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/lambda_support.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/bool_trait_undef.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/add_reference.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_reference.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/type_trait_def.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/type_trait_undef.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/ice.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/ice_or.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/ice_and.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/ice_not.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/ice_eq.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_arithmetic.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_integral.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_float.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_abstract.hpp
PolyhedralSurf.o: /usr/include/boost/static_assert.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_class.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_union.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/remove_cv.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/broken_compiler_spec.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/cv_traits_impl.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/intrinsics.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/identity.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/na_spec.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/lambda_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/void_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/na.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/na_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/lambda_arity_param.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/arity.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/dtp.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/preprocessor/enum.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/preprocessor/def_params_tail.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/limits/arity.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/logical/and.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/logical/bitand.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/identity.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/facilities/identity.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/empty.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/arithmetic/add.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/arithmetic/dec.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/control/while.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/list/fold_left.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/list/detail/fold_left.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/control/expr_iif.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/list/adt.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/detail/is_binary.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/detail/check.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/logical/compl.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/list/fold_right.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/list/detail/fold_right.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/list/reverse.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/control/detail/while.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/tuple/elem.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/arithmetic/sub.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/eti.hpp
PolyhedralSurf.o: /usr/include/boost/concept_archetype.hpp
PolyhedralSurf.o: /usr/include/stdlib.h /usr/include/sys/types.h
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
