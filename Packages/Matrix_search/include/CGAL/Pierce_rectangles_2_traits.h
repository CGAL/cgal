#line 70 "pierce_traits.awi"
#line 18 "code_formatting.awi"
// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// 2-4-Piercing Axis-Parallel 2D-Rectangles
// ============================================================================

#line 74 "pierce_traits.awi"
#line 54 "code_formatting.awi"
#if ! (CGAL_PIERCE_RECTANGLES_2_TRAITS_H)
#define CGAL_PIERCE_RECTANGLES_2_TRAITS_H 1

#line 25 "pierce_traits.awi"
#include <CGAL/Point_2.h>

#line 46 "code_formatting.awi"
CGAL_BEGIN_NAMESPACE
#line 28 "pierce_traits.awi"

template < class _R >
struct Piercing_traits_cartesian {
  // types
  typedef _R                           R;
  typedef typename R::FT               FT;
  typedef Point_2< R >                 Point_2;
  // constructions
  typedef Infinity_distance_2< R >     Infinity_distance_2;

  // get object methods:
  Infinity_distance_2
  get_infinity_distance_2() const
  { return Infinity_distance_2(); }

}; // Piercing_traits_cartesian
#line 66 "pierce_traits.awi"
#line 50 "code_formatting.awi"
CGAL_END_NAMESPACE
#line 67 "pierce_traits.awi"


#endif // ! (CGAL_PIERCE_RECTANGLES_2_TRAITS_H)

#line 12 "code_formatting.awi"
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

