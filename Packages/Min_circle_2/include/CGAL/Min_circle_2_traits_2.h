// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Min_circle_2_traits_2.h
// package       : $CGAL_Package: Min_circle_2 $
// chapter       : Geometric Optimisation
//
// source        : web/Min_circle_2.aw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>, Bernd Gärtner
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: default traits class for 2D Smallest Enclosing Circle
// ============================================================================

#ifndef CGAL_MIN_CIRCLE_2_TRAITS_2_H
#define CGAL_MIN_CIRCLE_2_TRAITS_2_H

// includes
#ifndef CGAL_POINT_2_H
#  include <CGAL/Point_2.h>
#endif
#ifndef CGAL_OPTIMISATION_CIRCLE_2_H
#  include <CGAL/Optimisation_circle_2.h>
#endif
#ifndef CGAL_PREDICATES_ON_POINTS_2_H
#  include <CGAL/predicates_on_points_2.h>
#endif

CGAL_BEGIN_NAMESPACE

// Class declarations
// ==================
template < class Traits_ >
class Min_circle_2;

template < class K_ >
class Min_circle_2_traits_2;

// Class interface and implementation
// ==================================
template < class K_ >
class Min_circle_2_traits_2 {
  public:
    // types
    typedef  K_                              K;
    typedef  CGAL::Point_2<K>                Point;
    typedef  CGAL::Optimisation_circle_2<K>  Circle;

private:
    // data members
    Circle  circle;                                 // current circle

    // friends
    friend  class CGAL::Min_circle_2< CGAL::Min_circle_2_traits_2<K> >;

  public:
    // creation (use default implementations)
    // CGAL::Min_circle_2_traits_2( );
    // CGAL::Min_circle_2_traits_2( CGAL::Min_circle_2_traits_2<K> const&);

    // operations
    inline
    CGAL::Orientation
    orientation( const Point& p, const Point& q, const Point& r) const
    {
        return( CGAL::orientation( p, q, r));
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_MIN_CIRCLE_2_TRAITS_2_H

// ===== EOF ==================================================================
