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
// release_date  : 2000, December 10
// 
// source        : Memory.fw
// file          : memory.h
// package       : Kernel_basic (3.17)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 3.17
// revision_date : 10 Dec 2000 
// author(s)     : Michael Seel
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================


#ifndef CGAL_MEMORY_H
#define CGAL_MEMORY_H
#include <memory>

#ifdef CGAL_USE_LEDA
#include <LEDA/allocator.h>
#define CGAL_ALLOCATOR(t) leda_allocator< t >
#define CGAL_ALLOC leda_allocator
#else
#define CGAL_ALLOCATOR(t) std::allocator< t >
#define CGAL_ALLOC std::allocator
#endif

#endif // CGAL_MEMORY_H
