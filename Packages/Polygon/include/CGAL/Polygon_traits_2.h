// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-0.9-I-06 $
// release_date  : $CGAL_Date: 1998/03/11 $
//
// file          : include/CGAL/Polygon_traits_2.h
// source        :
// revision      : 1.8a
// revision_date : 13 Mar 1998
// author(s)     : Wieger Wesselink <wieger@cs.ruu.nl>
//
// coordinator   : Utrecht University
//
// ============================================================================

#ifndef CGAL_POLYGON_TRAITS_2_H
#define CGAL_POLYGON_TRAITS_2_H

#ifndef CGAL_DIRECTION_2_H
#include <CGAL/Direction_2.h>
#endif // CGAL_DIRECTION_2_H
#ifndef CGAL_ISO_RECTANGLE_2_H
#include <CGAL/Iso_rectangle_2.h>
#endif // CGAL_ISO_RECTANGLE_2_H
#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_POINT_2_H
#ifndef CGAL_SEGMENT_2_H
#include <CGAL/Segment_2.h>
#endif // CGAL_SEGMENT_2_H
#ifndef CGAL_TRIANGLE_2_H
#include <CGAL/Triangle_2.h>
#endif // CGAL_TRIANGLE_2_H
#ifndef CGAL_VECTOR_2_H
#include <CGAL/Vector_2.h>
#endif // CGAL_VECTOR_2_H

#ifndef CGAL_NUMBER_UTILS_H
#include <CGAL/number_utils.h>
#endif // CGAL_NUMBER_UTILS_H
#ifndef CGAL_PREDICATES_ON_POINTS_2_H
#include <CGAL/predicates_on_points_2.h>
#endif // CGAL_PREDICATES_ON_POINTS_2_H
#ifndef CGAL_PREDICATE_OBJECTS_ON_POINTS_2_H
#include <CGAL/predicate_objects_on_points_2.h>
#endif // CGAL_PREDICATE_OBJECTS_ON_POINTS_2_H
#ifndef CGAL_SEGMENT_2_SEGMENT_2_INTERSECTION_H
#include <CGAL/Segment_2_Segment_2_intersection.h>
#endif // CGAL_SEGMENT_2_SEGMENT_2_INTERSECTION_H

//-----------------------------------------------------------------------//
//                          CGAL_Polygon_traits_2_aux
//-----------------------------------------------------------------------//
// The class CGAL_Polygon_traits_2_aux is needed as a workaround for compilers
// that have problems with member functions with return type R::FT.
// The template parameter _Point is added to simplify the use of points
// with additional information. It is assumed that _Point inherits from
// CGAL_Point_2<_R>.

template <class _R, class _FT, class _Point>
class CGAL_Polygon_traits_2_aux : public _R
{
  public:
    typedef _FT                           FT;
    typedef _Point                        Point_2;
    typedef CGAL_Segment_2<_R>            Segment_2;
    typedef CGAL_Vector_2<_R>             Vector_2;

    typedef CGAL_p_Less_xy<Point_2>       Less_xy;
    typedef CGAL_p_Less_yx<Point_2>       Less_yx;

    bool lexicographically_xy_smaller(const Point_2& p, const Point_2& q) const
    {
      return CGAL_lexicographically_xy_smaller(p,q);
    }

    bool lexicographically_yx_smaller_or_equal(const Point_2& p,
                                               const Point_2& q) const
    {
      return CGAL_lexicographically_yx_smaller_or_equal(p,q);
    }

    FT cross_product_2(const Vector_2& p, const Vector_2& q) const
    {
      return p.x() * q.y() - q.x() * p.y();
      // there should be multiple versions of this function
      // to distinguish between cartesian and homogeneous coordinates
      // (for efficiency reasons)
    }

    FT determinant_2(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      return cross_product_2(p-q, p-r);
    }

    int sign(const FT& x) const
    {
      return CGAL_sign(x);
    }

    bool is_negative(const FT& x) const
    {
      return CGAL_is_negative(x);
    }

    bool do_intersect(const Point_2& p1,
                      const Point_2& q1,
                      const Point_2& p2,
                      const Point_2& q2) const
    {
      return CGAL_do_intersect(Segment_2(p1,q1), Segment_2(p2,q2));
    }

    CGAL_Orientation orientation(const Point_2& p,
                                 const Point_2& q,
                                 const Point_2& r) const
    {
       return CGAL_orientation(p, q, r);
    }

    CGAL_Comparison_result compare_x(const Point_2 &p, const Point_2 &q) const
    {
       return CGAL_compare_x(p, q);
    }

    CGAL_Comparison_result compare_y(const Point_2 &p, const Point_2 &q) const
    {
       return CGAL_compare_y(p, q);
    }

    bool have_equal_direction(const Vector_2& v1,
                              const Vector_2& v2 ) const
    {
       return CGAL_Direction_2<_R>(v1) == CGAL_Direction_2<_R>(v2);
    }
};

//-----------------------------------------------------------------------//
//                          CGAL_Polygon_traits_2
//-----------------------------------------------------------------------//
// The class CGAL_Polygon_traits_2 is a traits class for CGAL_Polygon_2.

template <class _R>
class CGAL_Polygon_traits_2 :
  public CGAL_Polygon_traits_2_aux<_R, typename _R::FT, CGAL_Point_2<_R> >
{
};

#endif // CGAL_POLYGON_TRAITS_2_H

