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

#ifndef __SGI_STL_PTHREAD_ALLOC_H
#define __SGI_STL_PTHREAD_ALLOC_H

#include <stl_config.h>
#if defined  (__STL_DEBUG) || defined (__STL_ASSERTIONS)
# include <stldebug.h>
#endif

#include <pthread_alloc>

#ifdef __STL_USE_NAMESPACES

using __STLPORT_STD::_Pthread_alloc_template;
using __STLPORT_STD::pthread_alloc;

#endif /* __STL_USE_NAMESPACES */


#endif /* __SGI_STL_PTHREAD_ALLOC_H */

// Local Variables:
// mode:C++
// End:
