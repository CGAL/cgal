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

#ifndef __SGI_STL_ALGOBASE_H
#define __SGI_STL_ALGOBASE_H

#ifndef __SGI_STL_PAIR_H
#include <pair.h>
#endif

// memmove
#ifndef __STLPORT_CSTRING
# include <cstring>
#endif

// CHAR_MAX
#ifndef __STLPORT_CLIMITS
# include <climits>
#endif

#ifndef __SGI_STL_ITERATOR_H
#include <iterator.h>
#endif

#ifndef __SGI_STL_INTERNAL_ALGOBASE_H
#include <stl_algobase.h>
#endif

#ifndef __SGI_STL_INTERNAL_UNINITIALIZED_H
#include <stl_uninitialized.h>
#endif

#ifdef __STL_USE_NAMESPACES

# ifdef __STL_BROKEN_USING_DIRECTIVE
using namespace __STLPORT_STD;
# else
// Names from stl_algobase.h
using __STLPORT_STD::iter_swap; 
using __STLPORT_STD::swap; 
using __STLPORT_STD::min; 
using __STLPORT_STD::max; 
using __STLPORT_STD::copy; 
using __STLPORT_STD::copy_backward; 
using __STLPORT_STD::copy_n; 
using __STLPORT_STD::fill; 
using __STLPORT_STD::fill_n; 
using __STLPORT_STD::mismatch; 
using __STLPORT_STD::equal; 
using __STLPORT_STD::lexicographical_compare; 
using __STLPORT_STD::lexicographical_compare_3way; 

// Names from stl_uninitialized.h
using __STLPORT_STD::uninitialized_copy;
using __STLPORT_STD::uninitialized_copy_n;
using __STLPORT_STD::uninitialized_fill;
using __STLPORT_STD::uninitialized_fill_n;
# endif /* __STL_BROKEN_USING_DIRECTIVE */
#endif /* __STL_USE_NAMESPACES */
#endif /* __SGI_STL_ALGOBASE_H */

// Local Variables:
// mode:C++
// End:
