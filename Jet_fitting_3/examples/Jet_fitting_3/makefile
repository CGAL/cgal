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
LAPACK_LDLIBS  = -L$(LAPACK_DIR) -L$(LAPACK_DIR)/F2CLIBS \
	-llapack_LINUX \
          -lblas_LINUX -ltmglib_LINUX \
          -lF77 -lI77 -lm -lc

#F2CDIR  = $(LAPACK_DIR)/F2CLIBS
#LAPACK_LDLIBS  = $(LAPACK_DIR)/lapack_LINUX.a \
#           $(LAPACK_DIR)/blas_LINUX.a \
#           $(F2CDIR)/libF77.a $(F2CDIR)/libI77.a -lm -lc

## it should not be needed $(LAPACK_DIR)/tmglib_LINUX.a \
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

PolyhedralSurf.o: /usr/include/boost/property_map.hpp
PolyhedralSurf.o: /usr/include/boost/config.hpp
PolyhedralSurf.o: /usr/include/boost/config/user.hpp
PolyhedralSurf.o: /usr/include/boost/config/select_compiler_config.hpp
PolyhedralSurf.o: /usr/include/boost/config/compiler/gcc.hpp
PolyhedralSurf.o: /usr/include/boost/config/select_stdlib_config.hpp
PolyhedralSurf.o: /usr/include/boost/config/select_platform_config.hpp
PolyhedralSurf.o: /usr/include/boost/config/suffix.hpp /usr/include/limits.h
PolyhedralSurf.o: /usr/include/features.h /usr/include/sys/cdefs.h
PolyhedralSurf.o: /usr/include/gnu/stubs.h
PolyhedralSurf.o: /usr/lib/gcc/i386-redhat-linux/3.4.4/include/limits.h
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
PolyhedralSurf.o: /usr/include/boost/mpl/bool.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/bool_fwd.hpp
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
PolyhedralSurf.o: /usr/include/boost/static_assert.hpp
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
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/cv_traits_impl.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_volatile.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/remove_const.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/broken_compiler_spec.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/remove_volatile.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/remove_cv.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/function_traits.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_function.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/false_result.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/detail/is_function_ptr_helper.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/add_pointer.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/remove_reference.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/iterator_categories.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/detail/config_def.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/eval_if.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/if.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/value_wknd.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/aux_/config/integral.hpp
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
PolyhedralSurf.o: /usr/include/boost/type_traits/is_class.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_union.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/intrinsics.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/remove_pointer.hpp
PolyhedralSurf.o: /usr/include/boost/mpl/not.hpp
PolyhedralSurf.o: /usr/include/boost/iterator/detail/enable_if.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_pod.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_void.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_scalar.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_arithmetic.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_integral.hpp
PolyhedralSurf.o: /usr/include/boost/type_traits/is_float.hpp
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
PolyhedralSurf.o: /usr/include/stdlib.h
PolyhedralSurf.o: /usr/lib/gcc/i386-redhat-linux/3.4.4/include/stddef.h
PolyhedralSurf.o: /usr/include/stdio.h /usr/include/bits/types.h
PolyhedralSurf.o: /usr/include/bits/wordsize.h /usr/include/bits/typesizes.h
PolyhedralSurf.o: /usr/include/libio.h /usr/include/_G_config.h
PolyhedralSurf.o: /usr/include/wchar.h /usr/include/bits/wchar.h
PolyhedralSurf.o: /usr/include/gconv.h
PolyhedralSurf.o: /usr/lib/gcc/i386-redhat-linux/3.4.4/include/stdarg.h
PolyhedralSurf.o: /usr/include/bits/stdio_lim.h
PolyhedralSurf.o: /usr/include/bits/sys_errlist.h
