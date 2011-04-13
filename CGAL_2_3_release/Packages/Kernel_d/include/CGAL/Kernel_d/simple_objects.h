// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Kernel_d/simple_objects.h
// package       : Kernel_d
// chapter       : Basic
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Susan Hert <hert@mpi-sb.mpg.de>
//
// implementation: simple functor objects for the kernel
// ============================================================================
#ifndef SIMPLE_OBJECTS_H
#define SIMPLE_OBJECTS_H

CGAL_BEGIN_NAMESPACE

template <class R>
struct Lt_from_compare {
  typedef typename R::Point_d Point_d;
  bool operator()(const Point_d& p1, const Point_d& p2) const
  { typename R::Compare_lexicographically_d cmp;
    return cmp(p1,p2) == SMALLER; }
};

template <class R>
struct Le_from_compare {
  typedef typename R::Point_d Point_d;
  bool operator()(const Point_d& p1, const Point_d& p2) const
  { typename R::Compare_lexicographically_d cmp;
    return cmp(p1,p2) != LARGER; }
};

template <class R>
struct Eq_from_method {
  typedef typename R::Point_d Point_d;
  bool operator()(const Point_d& p1, const Point_d& p2) const
  { return p1 == p2; }
};

CGAL_END_NAMESPACE
#endif //SIMPLE_OBJECTS_H
