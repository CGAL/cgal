// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
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
// ======================================================================

#ifndef CGAL_POLYGON_TRAITS_2_H
#define CGAL_POLYGON_TRAITS_2_H

/*
#include <CGAL/Direction_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Vector_2.h>

#include <CGAL/number_utils.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/predicate_classes_2.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>
*/

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------//
//                          Polygon_traits_2_aux
//-----------------------------------------------------------------------//
// The class Polygon_traits_2_aux is needed as a workaround for compilers
// that have problems with member functions with return type R::FT.
// The template parameter _Point is added to simplify the use of points
// with additional information. It is assumed that _Point inherits from
// Point_2<_R>.
/*
template <class _R, class _FT, class _Point>
class Polygon_traits_2_aux : public _R
{
  public:
    typedef _FT                           FT;
    typedef _Point                        Point_2;
    typedef ::CGAL::Segment_2<_R>            Segment_2;
    typedef ::CGAL::Vector_2<_R>             Vector_2;

    typedef ::CGAL::Less_xy_2<Point_2>       Less_xy;
    typedef ::CGAL::Less_yx_2<Point_2>       Less_yx;

    bool lexicographically_xy_smaller(const Point_2& p, const Point_2& q) const
    {
      return ::CGAL::lexicographically_xy_smaller(p,q);
    }

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

    FT determinant_2(const Point_2& p,
                     const Point_2& q, const Point_2& r) const
    {
      return cross_product_2(p-q, p-r);
    }

    int sign(const FT& x) const
    {
      return CGAL_NTS sign(x);
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

    Orientation orientation(const Point_2& p,
                                 const Point_2& q,
                                 const Point_2& r) const
    {
       return ::CGAL::orientation(p, q, r);
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
       return ::CGAL::Direction_2<_R>(v1) == ::CGAL::Direction_2<_R>(v2);
    }
};
*/
//-----------------------------------------------------------------------//
//                          Polygon_traits_2
//-----------------------------------------------------------------------//
// The class Polygon_traits_2 is a traits class for Polygon_2.
/*
template <class _R>
class Polygon_traits_2 :
  public Polygon_traits_2_aux<_R, typename _R::FT, Point_2<_R> >
{
};
*/
template <class R_>
class Polygon_traits_2 : public R_ {};

CGAL_END_NAMESPACE

#endif // CGAL_POLYGON_TRAITS_2_H

