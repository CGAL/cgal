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
// file          : include/CGAL/Optimisation/Construct_point_d.h
// package       : $CGAL_Package: Optimisation_basic $
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: dD construction `point'
// ============================================================================

#ifndef CGAL_OPTIMISATION_CONSTRUCT_POINT_D_H
#define CGAL_OPTIMISATION_CONSTRUCT_POINT_D_H

// includes
#ifndef CGAL_POINT_D_H
#  include <CGAL/Point_d.h>
#endif

CGAL_BEGIN_NAMESPACE

// Class declaration
// =================
template < class R >
class Construct_point_d;

// Class interface
// ===============
template < class R_ >
class Construct_point_d {
  public:
    // self
    typedef  R_                         R;
    typedef  Construct_point_d<R>       Self;

    // types
    typedef  typename R::Point_d        Point;

    // creation
    Construct_point_d( ) { }

    // operations
    template < class InputIterator >
    Point
    operator() ( int d, InputIterator first, InputIterator last) const
    {
	return Point( d, first, last);
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_OPTIMISATION_CONSTRUCT_POINT_D_H

// ===== EOF ==================================================================
