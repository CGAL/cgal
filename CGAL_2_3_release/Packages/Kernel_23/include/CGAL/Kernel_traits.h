// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : include/CGAL/Kernel_traits.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_KERNEL_TRAITS_H
#define CGAL_KERNEL_TRAITS_H

CGAL_BEGIN_NAMESPACE

template <class T>
struct Kernel_traits
{
  typedef typename T::R Kernel;
};

CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_TRAITS_H
