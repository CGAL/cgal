// ============================================================================
//
// Copyright (c) 1998, 1999, 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : Pierce_rectangles_2_traits.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : pcenter.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// 2-4-Piercing Axis-Parallel 2D-Rectangles
// ============================================================================

#if ! (CGAL_PIERCE_RECTANGLES_2_TRAITS_H)
#define CGAL_PIERCE_RECTANGLES_2_TRAITS_H 1

#include <CGAL/Point_2.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
struct Piercing_traits_cartesian {
  // types
  typedef R_                           R;
  typedef typename R::FT               FT;
  typedef Point_2< R >                 Point_2;
  // constructions
  typedef Infinity_distance_2< R >     Infinity_distance_2;

  // get object methods:
  Infinity_distance_2
  get_infinity_distance_2() const
  { return Infinity_distance_2(); }

}; // Piercing_traits_cartesian
CGAL_END_NAMESPACE


#endif // ! (CGAL_PIERCE_RECTANGLES_2_TRAITS_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

