/*
 * Copyright (c) 1996,1997
 * Silicon Graphics Computer Systems, Inc.
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

#ifndef __SGI_STL_ALLOC_H
#define __SGI_STL_ALLOC_H

#ifndef __STL_CONFIG_H
#include <stl_config.h>
#endif

#if defined  (__STL_DEBUG) || defined (__STL_ASSERTIONS) && !defined (__STLPORT_DEBUG_H)
# include <stldebug.h>
#endif

# ifndef __STLPORT_CSTDDEF
#  include <cstddef>
# endif
# ifndef __STLPORT_CLIMITS
#  include <climits>
# endif
# ifndef __STLPORT_CSTDLIB
#  include <cstdlib>
# endif
# ifndef __STLPORT_CSTRING
#  include <cstring>
# endif
# ifndef __STLPORT_CASSERT
#  include <cassert>
# endif

#ifndef __SGI_STL_INTERNAL_ALLOC_H
#include <stl_alloc.h>
#endif

// Old SGI names
__STL_BEGIN_NAMESPACE

typedef __sgi_alloc alloc;
typedef __malloc_alloc<0> malloc_alloc;
#ifdef __STL_USE_NEWALLOC
typedef __new_alloc new_alloc;
#endif

#define simple_alloc __simple_alloc
typedef __single_client_alloc  single_client_alloc; 
typedef __multithreaded_alloc  multithreaded_alloc; 

__STL_END_NAMESPACE

#ifdef __STL_USE_NAMESPACES
# ifdef __STL_BROKEN_USING_DIRECTIVE

using namespace __STLPORT_STD;

# else

using __STLPORT_STD::malloc_alloc; 
using __STLPORT_STD::simple_alloc; 
# ifdef __STL_DEBUG_ALLOC
using __STLPORT_STD::__debug_alloc;
# endif 
#ifdef __STL_USE_NEWALLOC
using __STLPORT_STD::new_alloc;
#endif
using __STLPORT_STD::alloc; 
using __STLPORT_STD::single_client_alloc; 
using __STLPORT_STD::multithreaded_alloc; 
using __STLPORT_STD::allocator;

# endif /* __STL_BROKEN_USING_DIRECTIVE */
#endif /* __STL_USE_NAMESPACES */

#endif /* __SGI_STL_ALLOC_H */

// Local Variables:
// mode:C++
// End:
