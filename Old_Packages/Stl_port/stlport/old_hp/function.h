/*
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Copyright (c) 1996,1997
 * Silicon Graphics Computer Systems, Inc.
 *
 * Copyright (c) 1997
 * Moscow Center for SPARC Technology
 *
 * Copyright (c) 1999 
 * Boris Fomitchev
 *
 * This material is provided "as is", with absolutely no warranty expressed
 * or implied. Any use is at your own risk.
 *
 * Permission to use or copy this software for any purpose is hereby granted 
 * without fee, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is granted,
 * provided the above notices are retained, and a notice that the code was
 * modified is included with the above copyright notice.
 *
 */

#ifndef __SGI_STL_FUNCTION_H
#define __SGI_STL_FUNCTION_H

#ifndef __STL_CONFIG_H
#include <stl_config.h>
#endif

#ifndef __STLPORT_CSTDDEF
# include <cstddef>
#endif

#ifndef __SGI_STL_INTERNAL_RELOPS
#include <stl_relops.h>
#endif

#ifndef __SGI_STL_INTERNAL_FUNCTION_H
#include <stl_function.h>
#endif

#ifdef __STL_USE_NAMESPACES

# ifdef __STL_BROKEN_USING_DIRECTIVE
using namespace __STLPORT_STD;
#ifdef __STL_USE_NAMESPACE_FOR_RELOPS
using namespace __STLPORT_STD_RELOPS;
#endif /* __STL_USE_NAMESPACE_FOR_RELOPS */

# else /* __STL_BROKEN_USING_DIRECTIVE */

// Names from stl_function.h
using __STLPORT_STD::unary_function; 
using __STLPORT_STD::binary_function; 
using __STLPORT_STD::plus; 
using __STLPORT_STD::minus; 
using __STLPORT_STD::multiplies; 
using __STLPORT_STD::divides; 
using __STLPORT_STD::identity_element; 
using __STLPORT_STD::modulus; 
using __STLPORT_STD::negate; 
using __STLPORT_STD::equal_to; 
using __STLPORT_STD::not_equal_to; 
using __STLPORT_STD::greater; 
using __STLPORT_STD::less; 
using __STLPORT_STD::greater_equal; 
using __STLPORT_STD::less_equal; 
using __STLPORT_STD::logical_and; 
using __STLPORT_STD::logical_or; 
using __STLPORT_STD::logical_not; 
using __STLPORT_STD::unary_negate; 
using __STLPORT_STD::binary_negate; 
using __STLPORT_STD::not1; 
using __STLPORT_STD::not2; 
using __STLPORT_STD::binder1st; 
using __STLPORT_STD::binder2nd; 
using __STLPORT_STD::bind1st; 
using __STLPORT_STD::bind2nd; 
using __STLPORT_STD::unary_compose; 
using __STLPORT_STD::binary_compose; 
using __STLPORT_STD::compose1; 
using __STLPORT_STD::compose2; 
using __STLPORT_STD::pointer_to_unary_function; 
using __STLPORT_STD::pointer_to_binary_function; 
using __STLPORT_STD::ptr_fun; 
using __STLPORT_STD::identity; 
using __STLPORT_STD::select1st; 
using __STLPORT_STD::select2nd; 
using __STLPORT_STD::project1st; 
using __STLPORT_STD::project2nd; 
using __STLPORT_STD::constant_void_fun; 
using __STLPORT_STD::constant_unary_fun; 
using __STLPORT_STD::constant_binary_fun; 
using __STLPORT_STD::constant0; 
using __STLPORT_STD::constant1; 
using __STLPORT_STD::constant2; 
using __STLPORT_STD::subtractive_rng; 
using __STLPORT_STD::mem_fun_t; 
using __STLPORT_STD::const_mem_fun_t; 
using __STLPORT_STD::mem_fun_ref_t; 
using __STLPORT_STD::const_mem_fun_ref_t; 
using __STLPORT_STD::mem_fun1_t; 
using __STLPORT_STD::const_mem_fun1_t; 
using __STLPORT_STD::mem_fun1_ref_t; 
using __STLPORT_STD::const_mem_fun1_ref_t; 
using __STLPORT_STD::mem_fun; 
using __STLPORT_STD::mem_fun_ref; 
using __STLPORT_STD::mem_fun1; 
using __STLPORT_STD::mem_fun1_ref; 
# endif
#endif /* __STL_USE_NAMESPACES */

#endif /* __SGI_STL_FUNCTION_H */

// Local Variables:
// mode:C++
// End:
