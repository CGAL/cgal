// ======================================================================
//
// Copyright (c) 1999, 2002 The CGAL Consortium
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

// CGAL_ALLOCATOR(t) defines the default allocator used by CGAL.
// CGAL_MEMORY(t) overloads the new and delete operators for a given class.

#ifdef CGAL_USE_LEDA
#  include <LEDA/allocator.h>
#  ifndef CGAL_ALLOCATOR
#    define CGAL_ALLOCATOR(t) leda_allocator< t >
#  endif
#  ifndef CGAL_MEMORY
#    define CGAL_MEMORY(t)    LEDA_MEMORY(t)
#  endif
#else
#  include <memory>
#  ifndef CGAL_ALLOCATOR
#    define CGAL_ALLOCATOR(t) std::allocator< t >
// For debugging with GCC 3.1, the following allocator can be useful :
// std::__allocator<t, std::__debug_alloc<std::__malloc_alloc_template<0> > >
#  endif
#  ifndef CGAL_MEMORY
#    define CGAL_MEMORY(t)
#  endif
#endif

#endif // CGAL_MEMORY_H
