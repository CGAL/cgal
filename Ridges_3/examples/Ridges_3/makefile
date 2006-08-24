# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
include $(CGAL_MAKEFILE)

PROFOPT=  -O3 #-g 
CFLAGS=${PROFOPT} 
#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
LAPACK_INCDIRS = -I$(LAPACK_DIR)/SRC -I$(LAPACK_DIR)

#note line -I../../../Jet_fitting_3/include needed since 
#jet_fitting_3 not in the release yet
CXXFLAGS = \
           -I../../include \
           -I../../../Jet_fitting_3/include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   $(LAPACK_INCDIRS) 

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#
F2CDIR  = $(LAPACK_DIR)/F2CLIBS
LAPACK_LDLIBS  = $(LAPACK_DIR)/lapack_LINUX.a \
          $(LAPACK_DIR)/blas_LINUX.a $(LAPACK_DIR)/tmglib_LINUX.a \
          $(F2CDIR)/libF77.a $(F2CDIR)/libI77.a -lm -lc

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
	${LAPACK_LDLIBS} \
	$(LONG_NAME_PROBLEM_LDFLAGS) \
        $(CGAL_LDFLAGS) 

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                blind 

BLIND_OBJS=options.o PolyhedralSurf.o blind.o
blind: $(BLIND_OBJS)
	$(CGAL_CXX) $(LIBPATH) -o blind $(BLIND_OBJS) \
	$(LDFLAGS) 

clean: \
        blind.clean  options.clean PolyhedralSurf.clean 

depend: 
	makedepend *.[Ch]

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

# DO NOT DELETE

blind.o: /usr/include/boost/property_map.hpp /usr/include/boost/config.hpp
blind.o: /usr/include/boost/config/user.hpp
blind.o: /usr/include/boost/config/select_compiler_config.hpp
blind.o: /usr/include/boost/config/compiler/gcc.hpp
blind.o: /usr/include/boost/config/select_stdlib_config.hpp
blind.o: /usr/include/boost/config/select_platform_config.hpp
blind.o: /usr/include/boost/config/suffix.hpp /usr/include/limits.h
blind.o: /usr/include/features.h /usr/include/sys/cdefs.h
blind.o: /usr/include/gnu/stubs.h /usr/include/bits/posix1_lim.h
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
blind.o: ../../../Jet_fitting_3/include/CGAL/Monge_via_jet_fitting.h
blind.o: /usr/include/math.h /usr/include/bits/huge_val.h
blind.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
blind.o: ../../../Jet_fitting_3/include/CGAL/Lapack/Linear_algebra_lapack.h
blind.o: /usr/include/stdlib.h /usr/include/sys/types.h
blind.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
blind.o: /usr/include/bits/typesizes.h /usr/include/time.h
blind.o: /usr/include/endian.h /usr/include/bits/endian.h
blind.o: /usr/include/sys/select.h /usr/include/bits/select.h
blind.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
blind.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
blind.o: /usr/include/bits/sched.h /usr/include/alloca.h PolyhedralSurf.h
blind.o: /usr/include/stdio.h /usr/include/libio.h /usr/include/_G_config.h
blind.o: /usr/include/wchar.h /usr/include/bits/wchar.h /usr/include/gconv.h
blind.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
blind.o: PolyhedralSurf_operations.h PolyhedralSurf_rings.h options.h
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
PolyhedralSurf.o: /usr/include/bits/sys_errlist.h PolyhedralSurf_operations.h
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
PolyhedralSurf.o: /usr/include/bits/sys_errlist.h PolyhedralSurf_operations.h
