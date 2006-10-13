# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
include $(CGAL_MAKEFILE)

PROFOPT= -g #-O3 #-g
CFLAGS=${PROFOPT} 
#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
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
#this gives linking problems?? so i use the other
#LAPACK_LDLIBS  = -L$(LAPACK_DIR) -L$(LAPACK_DIR)/F2CLIBS \
#	-llapack_LINUX \
#          -lblas_LINUX -ltmglib_LINUX \
#          -lF77 -lI77 -lm -lc

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
                blind blind_1pt

BLIND_OBJS=options.o PolyhedralSurf.o blind.o
blind: $(BLIND_OBJS)
	$(CGAL_CXX) $(LIBPATH) -g -o blind $(BLIND_OBJS) \
	$(LDFLAGS) 

blind_1pt: blind_1pt.o
	$(CGAL_CXX) $(LIBPATH) -o blind_1pt blind_1pt.o \
	$(LDFLAGS) 

clean: \
        blind.clean blind_1pt.clean options.clean PolyhedralSurf.clean 

depend: 
	makedepend *.[Ch]

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
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
blind_1pt.o: /usr/include/alloca.h
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
blind.o: /usr/include/boost/concept_archetype.hpp PolyhedralSurf.h
blind.o: /usr/include/boost/graph/properties.hpp
blind.o: /usr/include/boost/pending/property.hpp
blind.o: /usr/include/boost/pending/ct_if.hpp
blind.o: /usr/include/boost/pending/detail/property.hpp
blind.o: /usr/include/boost/type_traits/same_traits.hpp
blind.o: /usr/include/boost/type_traits/is_same.hpp
blind.o: /usr/include/boost/graph/graph_traits.hpp
blind.o: /usr/include/boost/tuple/tuple.hpp /usr/include/boost/ref.hpp
blind.o: /usr/include/boost/utility/addressof.hpp
blind.o: /usr/include/boost/tuple/detail/tuple_basic.hpp
blind.o: /usr/include/boost/type_traits/cv_traits.hpp
blind.o: /usr/include/boost/type_traits/add_const.hpp
blind.o: /usr/include/boost/type_traits/add_volatile.hpp
blind.o: /usr/include/boost/type_traits/add_cv.hpp
blind.o: /usr/include/boost/type_traits/is_const.hpp
blind.o: /usr/include/boost/type_traits/is_volatile.hpp
blind.o: /usr/include/boost/type_traits/remove_const.hpp
blind.o: /usr/include/boost/type_traits/remove_volatile.hpp
blind.o: /usr/include/boost/type_traits/function_traits.hpp
blind.o: /usr/include/boost/type_traits/is_function.hpp
blind.o: /usr/include/boost/type_traits/detail/false_result.hpp
blind.o: /usr/include/boost/type_traits/detail/is_function_ptr_helper.hpp
blind.o: /usr/include/boost/type_traits/add_pointer.hpp
blind.o: /usr/include/boost/type_traits/remove_reference.hpp
blind.o: /usr/include/boost/mpl/if.hpp
blind.o: /usr/include/boost/mpl/aux_/value_wknd.hpp
blind.o: /usr/include/boost/mpl/aux_/config/integral.hpp
blind.o: /usr/include/boost/iterator/iterator_categories.hpp
blind.o: /usr/include/boost/iterator/detail/config_def.hpp
blind.o: /usr/include/boost/mpl/eval_if.hpp
blind.o: /usr/include/boost/mpl/placeholders.hpp
blind.o: /usr/include/boost/mpl/arg.hpp /usr/include/boost/mpl/arg_fwd.hpp
blind.o: /usr/include/boost/mpl/aux_/na_assert.hpp
blind.o: /usr/include/boost/mpl/aux_/arity_spec.hpp
blind.o: /usr/include/boost/mpl/aux_/arg_typedef.hpp
blind.o: /usr/include/boost/mpl/aux_/config/use_preprocessed.hpp
blind.o: /usr/include/boost/mpl/aux_/include_preprocessed.hpp
blind.o: /usr/include/boost/mpl/aux_/config/compiler.hpp
blind.o: /usr/include/boost/preprocessor/stringize.hpp
blind.o: /usr/include/boost/iterator/detail/config_undef.hpp
blind.o: /usr/include/boost/iterator/iterator_adaptor.hpp
blind.o: /usr/include/boost/iterator/iterator_facade.hpp
blind.o: /usr/include/boost/iterator/interoperable.hpp
blind.o: /usr/include/boost/mpl/or.hpp
blind.o: /usr/include/boost/mpl/aux_/nested_type_wknd.hpp
blind.o: /usr/include/boost/iterator/iterator_traits.hpp
blind.o: /usr/include/boost/iterator/detail/facade_iterator_category.hpp
blind.o: /usr/include/boost/mpl/and.hpp
blind.o: /usr/include/boost/detail/indirect_traits.hpp
blind.o: /usr/include/boost/type_traits/is_pointer.hpp
blind.o: /usr/include/boost/type_traits/is_member_pointer.hpp
blind.o: /usr/include/boost/type_traits/is_member_function_pointer.hpp
blind.o: /usr/include/boost/type_traits/detail/is_mem_fun_pointer_impl.hpp
blind.o: /usr/include/boost/type_traits/remove_pointer.hpp
blind.o: /usr/include/boost/mpl/not.hpp
blind.o: /usr/include/boost/iterator/detail/enable_if.hpp
blind.o: /usr/include/boost/implicit_cast.hpp
blind.o: /usr/include/boost/type_traits/is_pod.hpp
blind.o: /usr/include/boost/type_traits/is_void.hpp
blind.o: /usr/include/boost/type_traits/is_scalar.hpp
blind.o: /usr/include/boost/type_traits/is_enum.hpp
blind.o: /usr/include/boost/mpl/always.hpp /usr/include/boost/mpl/apply.hpp
blind.o: /usr/include/boost/mpl/apply_fwd.hpp
blind.o: /usr/include/boost/mpl/apply_wrap.hpp
blind.o: /usr/include/boost/mpl/aux_/has_apply.hpp
blind.o: /usr/include/boost/mpl/has_xxx.hpp
blind.o: /usr/include/boost/mpl/aux_/type_wrapper.hpp
blind.o: /usr/include/boost/mpl/aux_/yes_no.hpp
blind.o: /usr/include/boost/mpl/aux_/config/arrays.hpp
blind.o: /usr/include/boost/mpl/aux_/config/has_xxx.hpp
blind.o: /usr/include/boost/mpl/aux_/config/msvc_typename.hpp
blind.o: /usr/include/boost/mpl/aux_/config/has_apply.hpp
blind.o: /usr/include/boost/mpl/aux_/msvc_never_true.hpp
blind.o: /usr/include/boost/mpl/lambda.hpp /usr/include/boost/mpl/bind.hpp
blind.o: /usr/include/boost/mpl/bind_fwd.hpp
blind.o: /usr/include/boost/mpl/aux_/config/bind.hpp
blind.o: /usr/include/boost/mpl/next.hpp
blind.o: /usr/include/boost/mpl/next_prior.hpp
blind.o: /usr/include/boost/mpl/aux_/common_name_wknd.hpp
blind.o: /usr/include/boost/mpl/protect.hpp
blind.o: /usr/include/boost/mpl/aux_/full_lambda.hpp
blind.o: /usr/include/boost/mpl/quote.hpp /usr/include/boost/mpl/void.hpp
blind.o: /usr/include/boost/mpl/aux_/has_type.hpp
blind.o: /usr/include/boost/mpl/aux_/template_arity.hpp
blind.o: PolyhedralSurf_operations.h PolyhedralSurf_rings.h options.h
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
PolyhedralSurf.o: /usr/include/boost/graph/properties.hpp
PolyhedralSurf.o: /usr/include/boost/pending/property.hpp
PolyhedralSurf.o: /usr/include/boost/pending/ct_if.hpp
PolyhedralSurf.o: /usr/include/boost/pending/detail/property.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/same_traits.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_same.hpp
PolyhedralSurf.o: /usr/include/boost/graph/graph_traits.hpp
PolyhedralSurf.o: /usr/include/boost/tuple/tuple.hpp
PolyhedralSurf.o: /usr/include/boost/ref.hpp
PolyhedralSurf.o: /usr/include/boost/utility/addressof.hpp
PolyhedralSurf.o: /usr/include/boost/tuple/detail/tuple_basic.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/cv_traits.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/add_const.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/add_volatile.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/add_cv.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_const.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_volatile.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/remove_const.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/remove_volatile.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/function_traits.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_function.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/false_result.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/is_function_ptr_helper.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/add_pointer.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/remove_reference.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/if.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/value_wknd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/integral.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/iterator_categories.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/detail/config_def.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/eval_if.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/placeholders.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/arg.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/arg_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/na_assert.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/arity_spec.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/arg_typedef.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/use_preprocessed.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/include_preprocessed.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/compiler.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/stringize.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/detail/config_undef.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/iterator_adaptor.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/iterator_facade.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/interoperable.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/or.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/nested_type_wknd.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/iterator_traits.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/detail/facade_iterator_category.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/and.hpp
PolyhedralSurf.o: /usr/include/boost/detail/indirect_traits.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_pointer.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_member_pointer.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_member_function_pointer.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/is_mem_fun_pointer_impl.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/remove_pointer.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/not.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/detail/enable_if.hpp
PolyhedralSurf.o: /usr/include/boost/implicit_cast.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_pod.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_void.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_scalar.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_enum.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/always.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/apply.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/apply_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/apply_wrap.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/has_apply.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/has_xxx.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/type_wrapper.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/yes_no.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/arrays.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/has_xxx.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/msvc_typename.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/has_apply.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/msvc_never_true.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/lambda.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/bind.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/bind_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/bind.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/next.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/next_prior.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/common_name_wknd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/protect.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/full_lambda.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/quote.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/void.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/has_type.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/template_arity.hpp
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
PolyhedralSurf.o: /usr/include/boost/graph/properties.hpp
PolyhedralSurf.o: /usr/include/boost/pending/property.hpp
PolyhedralSurf.o: /usr/include/boost/pending/ct_if.hpp
PolyhedralSurf.o: /usr/include/boost/pending/detail/property.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/same_traits.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_same.hpp
PolyhedralSurf.o: /usr/include/boost/graph/graph_traits.hpp
PolyhedralSurf.o: /usr/include/boost/tuple/tuple.hpp
PolyhedralSurf.o: /usr/include/boost/ref.hpp
PolyhedralSurf.o: /usr/include/boost/utility/addressof.hpp
PolyhedralSurf.o: /usr/include/boost/tuple/detail/tuple_basic.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/cv_traits.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/add_const.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/add_volatile.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/add_cv.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_const.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_volatile.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/remove_const.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/remove_volatile.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/function_traits.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_function.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/false_result.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/is_function_ptr_helper.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/add_pointer.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/remove_reference.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/if.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/value_wknd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/integral.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/iterator_categories.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/detail/config_def.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/eval_if.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/placeholders.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/arg.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/arg_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/na_assert.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/arity_spec.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/arg_typedef.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/use_preprocessed.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/include_preprocessed.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/compiler.hpp
PolyhedralSurf.o: /usr/include/boost/preprocessor/stringize.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/detail/config_undef.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/iterator_adaptor.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/iterator_facade.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/interoperable.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/or.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/nested_type_wknd.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/iterator_traits.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/detail/facade_iterator_category.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/and.hpp
PolyhedralSurf.o: /usr/include/boost/detail/indirect_traits.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_pointer.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_member_pointer.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_member_function_pointer.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/is_mem_fun_pointer_impl.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/remove_pointer.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/not.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/detail/enable_if.hpp
PolyhedralSurf.o: /usr/include/boost/implicit_cast.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_pod.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_void.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_scalar.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_enum.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/always.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/apply.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/apply_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/apply_wrap.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/has_apply.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/has_xxx.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/type_wrapper.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/yes_no.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/arrays.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/has_xxx.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/msvc_typename.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/has_apply.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/msvc_never_true.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/lambda.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/bind.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/bind_fwd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/bind.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/next.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/next_prior.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/common_name_wknd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/protect.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/full_lambda.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/quote.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/void.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/has_type.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/template_arity.hpp
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
