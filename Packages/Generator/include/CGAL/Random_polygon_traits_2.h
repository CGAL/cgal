// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.2-I-12 $
// release_date  : $CGAL_Date: 2000/04/07 $
//
// file          : include/CGAL/Random_polygon_traits_2.h
// package       : $CGAL_Package $
// chapter       : Geometric Object Generators
//
// revision      : 1.0
// revision_date : 19 April 2000
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// implementation: Random Simple Polygon Traits
// ======================================================================

#ifndef CGAL_RANDOM_POLYGON_TRAITS_2_H
#define CGAL_RANDOM_POLYGON_TRAITS_2_H

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------//
//                          Random_polygon_traits_2
//-----------------------------------------------------------------------//

template <class R_>
class Random_polygon_traits_2 
{
  public:
    typedef R_                            R;
    typedef typename R::FT                FT;
    typedef typename R::Point_2           Point_2;
    typedef typename R::Less_xy_2         Less_xy_2;
    typedef typename R::Orientation_2     Orientation_2;

    Less_xy_2
    less_xy_2_object() const
    { return Less_xy_2(); }

    Orientation_2
    orientation_2_object() const
    { return Orientation_2(); }
};

CGAL_END_NAMESPACE

#endif // CGAL_RANDOM_POLYGON_TRAITS_2_H

