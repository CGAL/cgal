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

#ifndef __SGI_STL_HASH_SET_H
#define __SGI_STL_HASH_SET_H

#ifndef __SGI_STL_HASHTABLE_H
#include <hashtable.h>
#endif 

#include <stl_hash_set.h>

#ifdef __STL_USE_NAMESPACES
# ifdef __STL_BROKEN_USING_DIRECTIVE
using namespace __STLPORT_STD;
# else
using __STLPORT_STD::hash;
using __STLPORT_STD::hashtable;
using __STLPORT_STD::hash_set;
using __STLPORT_STD::hash_multiset;
using __STLPORT_STD::__hash_set__;
using __STLPORT_STD::__hash_multiset__;
# endif
#endif /* __STL_USE_NAMESPACES */

#endif /* __SGI_STL_HASH_SET_H */
