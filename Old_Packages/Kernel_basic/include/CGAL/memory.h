// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// author(s)     : Michael Seel
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_MEMORY_H
#define CGAL_MEMORY_H

// This file defines the macro CGAL_ALLOCATOR(t) which is the default
// allocator used by CGAL.

#ifdef CGAL_USE_LEDA
#  include <LEDA/allocator.h>
#  define CGAL_ALLOCATOR(t) leda_allocator< t >

#else
#  include <memory>
#  define CGAL_ALLOCATOR(t) std::allocator< t >
#endif

#endif // CGAL_MEMORY_H
