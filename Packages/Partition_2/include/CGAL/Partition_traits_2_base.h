// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : include/CGAL/Partition_traits_2_base.h
// package       : $CGAL_Package: Partition_2 $
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
// chapter       : Planar Polygon Partitioning
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: Base class for polygon partitioning function traits classes
// ============================================================================

#ifndef PARTITION_TRAITS_2_BASE_H
#define PARTITION_TRAITS_2_BASE_H

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_traits_2.h>
#include <list>

namespace CGAL {

template <class R_>
class Partition_traits_2_base
{
  public:
    typedef R_                                      R;
    typedef CGAL::Polygon_traits_2<R>               Poly_Traits;
    typedef typename Poly_Traits::Point_2           Point_2;
    typedef ::std::list<Point_2>                    Container;
    typedef CGAL::Polygon_2<Poly_Traits, Container> Polygon_2;
    typedef typename R::Less_yx_2                   Less_yx_2;
    typedef typename R::Less_xy_2                   Less_xy_2;
    typedef typename R::Leftturn_2                  Leftturn_2;
    typedef typename R::Orientation_2               Orientation_2;
    typedef typename R::Compare_y_2                 Compare_y_2;
    typedef typename R::Compare_x_2                 Compare_x_2;

    // needed below for do_intersect
    typedef typename R::Segment_2                   Segment_2;

    Less_yx_2
    less_yx_2_object() const
    { return Less_yx_2(); }

    Less_xy_2
    less_xy_2_object() const
    { return Less_xy_2(); }

    Leftturn_2
    leftturn_2_object() const
    { return Leftturn_2(); }

    Orientation_2
    orientation_2_object() const
    { return Orientation_2(); }

    Compare_y_2
    compare_y_2_object() const
    {  return Compare_y_2(); }

    Compare_x_2
    compare_x_2_object() const
    {  return Compare_x_2(); }


    // Needed for Polygon_algorithms_2 without the _2
    typedef typename R::Less_xy_2            Less_xy;

    // Needed for simplicity and CCW order precondition checking
    // ----------
    typedef typename Poly_Traits::Vector_2   Vector_2;
    typedef typename R::FT                   FT;

    typedef typename R::Direction_2          Direction_2;

    Comparison_result 
    compare_x(const Point_2 &p, const Point_2 &q) const
    {
       return Compare_x_2()(p, q);
    }

    Comparison_result 
    compare_y(const Point_2 &p, const Point_2 &q) const
    {
       return Compare_y_2()(p, q);
    }


    FT 
    cross_product_2(const Vector_2& p, const Vector_2& q) const
    {
        return p.x() * q.y() - q.x() * p.y();
    }

    bool 
    do_intersect(const Point_2& p1, const Point_2& q1, const Point_2& p2,
                 const Point_2& q2) const
    {
       return ::CGAL::do_intersect(Segment_2(p1,q1), Segment_2(p2,q2));
    }

    bool 
    have_equal_direction(const Vector_2& v1, const Vector_2& v2) const
    {
       return Direction_2(v1) == Direction_2(v2);
    }

    bool 
    is_negative(const FT& x) const
    {
       return CGAL_NTS is_negative(x);
    }

    // necessary for precondition checking
    bool 
    lexicographically_yx_smaller_or_equal(const Point_2& p, 
                                          const Point_2& q) const
    {
       return ::CGAL::lexicographically_yx_smaller_or_equal(p,q);
    }
    // ----------

    // necessary for postcondition and precondition checking
    ::CGAL::Orientation 
    orientation(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
       return Orientation_2()(p,q,r);
    }


};

}

#endif // PARTITION_TRAITS_2_BASE_H
