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
// file          : include/CGAL/Optimisation/Access_dimension_d.h
// package       : $CGAL_Package: Optimisation_basic $
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: dD data accessor `dimension'
// ============================================================================

#ifndef CGAL_OPTIMISATION_ACCESS_DIMENSION_D_H
#define CGAL_OPTIMISATION_ACCESS_DIMENSION_D_H

CGAL_BEGIN_NAMESPACE

// Class declaration
// =================
template < class R_ >
class Access_dimension_d;

// Class interface
// ===============
template < class R_ >
class Access_dimension_d {
  public:
    // self
    typedef  R_                         R;
    typedef  Access_dimension_d<R>      Self;

    // types
    typedef  typename R::Point_d        Point;

    // unary function class types
    typedef  int                        result_type;
    typedef  Point                      argument_type;

    // creation
    Access_dimension_d( ) { }

    // operations
    int  operator() ( const Point& p) const { return p.dimension(); }
};

CGAL_END_NAMESPACE

#endif // CGAL_OPTIMISATION_ACCESS_DIMENSION_D_H

// ===== EOF ==================================================================
