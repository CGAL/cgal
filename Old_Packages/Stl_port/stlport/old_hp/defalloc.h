/*
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Hewlett-Packard Company makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 */

// Inclusion of this file is DEPRECATED.  This is the original HP
// default allocator.  It is provided only for backward compatibility.
// This file WILL BE REMOVED in a future release.
//
// DO NOT USE THIS FILE unless you have an old container implementation
// that requires an allocator with the HP-style interface.  
//
// Standard-conforming allocators have a very different interface.  The
// standard default allocator is declared in the header <memory>.

// Adaptation note: THIS version of allocator<T> is fully compatible with
// SGI containers and works OK standalone. It is also as close to CD2 version
// as possible w/o member templates.
// However, explicit use of allocator<T>  is not recommended 
// unless you have to do so ( for example, compiling third-party code).

#ifndef __SGI_STL_DEFALLOC_H
#define __SGI_STL_DEFALLOC_H

# ifndef __SGI_STL_CONFIG_H
#  include <stl_config.h>
# endif

# if defined (__STL_DEBUG) && ! defined ( __SGI_STL_DEBUG_H )
#  include <stldebug.h>
# endif

#if defined (__STL_USE_NEW_STYLE_HEADERS)
# include <cstddef>
# include <cstdlib>
# include <cstring>
# include <cassert>
#else
# include <stddef.h>
# include <stdlib.h>
# include <string.h>
# include <assert.h>
#endif

# include <new>

#ifdef __STL_THREADS
# include <stl_threads.h>
#endif

# if !defined (__THROW_BAD_ALLOC) && !defined(__STL_USE_EXCEPTIONS)
#   if defined (__STL_USE_NEW_STYLE_HEADERS)
#    include <cstdio>
#   else
#    include <stdio.h>
#   endif
# endif
#  include <stl_alloc.h>

// fbp: just for backwards compatibility,
// hope this doesn't break anything.
#ifdef __STL_USE_NAMESPACES
# ifdef __STL_BROKEN_USING_DIRECTIVE
using namespace __STLPORT_STD;
# else
using __STLPORT_STD::allocator;
# endif /* __STL_BROKEN_USING_DIRECTIVE */
#endif /*  __STL_USE_NAMESPACES */
#endif /* __SGI_STL_DEFALLOC_H */
