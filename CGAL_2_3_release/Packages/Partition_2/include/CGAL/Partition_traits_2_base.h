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

};

}

#endif // PARTITION_TRAITS_2_BASE_H
