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

#ifndef __SGI_STL_STACK_H
#define __SGI_STL_STACK_H

#ifndef __SGI_STL_VECTOR_H
#include <vector.h>
#endif

#ifndef __SGI_STL_DEQUE_H
#include <deque.h>
#endif

#ifndef __SGI_STL_HEAP_H
#include <heap.h>
#endif

#ifndef __SGI_STL_INTERNAL_STACK_H
#include <stl_stack.h>
#endif

#ifndef __SGI_STL_INTERNAL_QUEUE_H
#include <stl_queue.h>
#endif

#ifdef __STL_USE_NAMESPACES
# ifdef __STL_BROKEN_USING_DIRECTIVE
using namespace __STLPORT_STD;
# else
using __STLPORT_STD::stack;
using __STLPORT_STD::queue;
using __STLPORT_STD::priority_queue;
# endif
#endif /* __STL_USE_NAMESPACES */

#endif /* __SGI_STL_STACK_H */

// Local Variables:
// mode:C++
// End:
