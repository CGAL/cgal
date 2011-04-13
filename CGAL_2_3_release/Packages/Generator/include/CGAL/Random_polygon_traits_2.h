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

#include <CGAL/Direction_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Vector_2.h>

#include <CGAL/number_utils.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/predicate_classes_2.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>

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
    typedef typename R::Segment_2         Segment_2;
    typedef typename R::Direction_2       Direction_2;
    typedef typename R::Vector_2          Vector_2;
    typedef typename R::Less_yx_2         Less_yx_2;
    typedef typename R::Less_xy_2         Less_xy_2;
    typedef typename R::Orientation_2     Orientation_2;

    Less_xy_2
    less_xy_2_object()
    { return Less_xy_2(); }

    Orientation_2
    orientation_2_object()
    { return Orientation_2(); }

    bool lexicographically_yx_smaller_or_equal(const Point_2& p,
                                               const Point_2& q) const
    {
      return ::CGAL::lexicographically_yx_smaller_or_equal(p,q);
    }

    FT cross_product_2(const Vector_2& p, const Vector_2& q) const
    {
      return p.x() * q.y() - q.x() * p.y();
      // there should be multiple versions of this function
      // to distinguish between cartesian and homogeneous coordinates
      // (for efficiency reasons)
    }


    bool is_negative(const FT& x) const
    {
      return CGAL_NTS is_negative(x);
    }

    bool do_intersect(const Point_2& p1,
                      const Point_2& q1,
                      const Point_2& p2,
                      const Point_2& q2) const
    {
      return ::CGAL::do_intersect(Segment_2(p1,q1), Segment_2(p2,q2));
    }

    Comparison_result compare_x(const Point_2 &p, const Point_2 &q) const
    {
       return ::CGAL::compare_x(p, q);
    }

    Comparison_result compare_y(const Point_2 &p, const Point_2 &q) const
    {
       return ::CGAL::compare_y(p, q);
    }

    bool have_equal_direction(const Vector_2& v1,
                              const Vector_2& v2 ) const
    {
       return Direction_2(v1) == Direction_2(v2);
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_RANDOM_POLYGON_TRAITS_2_H

