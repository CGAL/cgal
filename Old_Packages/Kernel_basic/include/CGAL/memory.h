// ======================================================================
//
// Copyright (c) 1999, 2002, 2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : memory.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Seel, Sylvain Pion
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_MEMORY_H
#define CGAL_MEMORY_H

#include <memory>

// CGAL_ALLOCATOR(t) defines the default allocator used by CGAL.
// CGAL_MEMORY(t) overloads the new and delete operators for a given class.

// When LEDA is there, the user could define these macros as
// leda_allocator< T >  and  LEDA_MEMORY(T) for example.

// For debugging with GCC, the following allocator can be useful :
// std::__allocator<T, std::__debug_alloc<std::__malloc_alloc_template<0> > >

#ifndef CGAL_ALLOCATOR
#  define CGAL_ALLOCATOR(T) std::allocator< T >
#endif

#ifndef CGAL_MEMORY
#  define CGAL_MEMORY(T)
#endif

#endif // CGAL_MEMORY_H
