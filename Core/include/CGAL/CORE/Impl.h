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

// macros for defining namespace
#define CORE_BEGIN_NAMESPACE    namespace CORE {
#define CORE_END_NAMESPACE      }

// condition preprocessor for inline function
#ifndef CORE_DISABLE_INLINE
  #define CORE_INLINE inline
#else
  #define CORE_INLINE
#endif

// Macros for memory pool
#ifdef CORE_DISABLE_MEMORY_POOL
  #define CORE_MEMORY(T)
#else
  #include <CGAL/CORE/MemoryPool.h>
  #define CORE_MEMORY(T)                                                 \
    void *operator new( size_t size)                                     \
    { return MemoryPool<T>::global_allocator().allocate(size); }         \
    void operator delete( void *p, size_t )                              \
    { MemoryPool<T>::global_allocator().free(p); }
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

