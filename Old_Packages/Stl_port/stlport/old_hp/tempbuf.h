/*
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
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

#ifndef __SGI_STL_TEMPBUF_H
#define __SGI_STL_TEMPBUF_H

# ifndef __STLPORT_CLIMITS
#  include <climits>
# endif

#ifndef __STLPORT_CSTDDEF
# include <stddef>
#endif

#ifndef __STLPORT_CSTDLIB
# include <cstdlib>
#endif

#if !defined (__STLPORT_NEW)
# include <new>
#endif

#ifndef __TYPE_TRAITS_H
#include <type_traits.h>
#endif

#ifndef __SGI_STL_PAIR_H
#include <pair.h>
#endif

#ifndef __SGI_STL_INTERNAL_TEMPBUF_H
#include <stl_tempbuf.h>
#endif

#ifdef __STL_USE_NAMESPACES

# ifdef __STL_BROKEN_USING_DIRECTIVE
using namespace __STLPORT_STD;
# else
using __STLPORT_STD::get_temporary_buffer;
using __STLPORT_STD::return_temporary_buffer;
using __STLPORT_STD::temporary_buffer;
# endif
#endif /* __STL_USE_NAMESPACES */

#endif /* __SGI_STL_TEMPBUF_H */

// Local Variables:
// mode:C++
// End:
