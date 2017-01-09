/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 * You can redistribute it and/or modify it under the terms of the GNU
 * Lesser General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * $URL$
 * $Id$
 ***************************************************************************/

#ifndef _CORE_IMPL_H_
#define _CORE_IMPL_H_

#include <CGAL/CORE/Config.h>

// The following lines only for MS Visual C++
#ifdef _MSC_VER
  #pragma warning(disable: 4291) // no matching operator delete found
  #pragma warning(disable: 4146) 
  #pragma warning(disable: 4267)
  #pragma warning(disable: 4244) 
#endif

// condition preprocessor for inline function
#ifndef CORE_DISABLE_INLINE
  #define CORE_INLINE inline
#else
  #define CORE_INLINE
#endif

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

// Macros for memory pool
#ifdef CORE_DISABLE_MEMORY_POOL
  #define CORE_NEW(T)
  #define CORE_DELETE(T)
  #define CORE_MEMORY_IMPL(T)
  #define CORE_MEMORY_IMPL_TEMPLATE_WITH_ONE_ARG(T)
#else
  #include <CGAL/CORE/MemoryPool.h>
  #define CORE_NEW(T)  void *operator new( size_t size);
  #define CORE_DELETE(T)  void operator delete( void *p, size_t );

  #define CORE_MEMORY_IMPL(T)                                            \
    CGAL_INLINE_FUNCTION void *T::operator new( size_t size)             \
    { return MemoryPool<T>::global_allocator().allocate(size); }         \
    CGAL_INLINE_FUNCTION void T::operator delete( void *p, size_t )      \
    { MemoryPool<T>::global_allocator().free(p); }
  #define CORE_MEMORY_IMPL_TEMPLATE_WITH_ONE_ARG(C)                       \
  template <typename T>                                                   \
  CGAL_INLINE_FUNCTION void *C<T>::operator new( size_t size)             \
  { return MemoryPool<C<T> >::global_allocator().allocate(size); }        \
  template <typename T>                                                   \
  CGAL_INLINE_FUNCTION void C<T>::operator delete( void *p, size_t )      \
  { MemoryPool<C<T> >::global_allocator().free(p); }
#endif

// include some common header files
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <cctype>
#include <climits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#endif // _CORE_IMPL_H_
