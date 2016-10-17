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
// release       : $CGAL_Revision: CGAL-2.5-I-119 $
// release_date  : $CGAL_Date: 2003/08/01 $
//
// file          : include/CGAL/Partition_traits_2.h
// package       : Partition_2 (1.53)
// maintainer    : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
// chapter       : Planar Polygon Partitioning
//
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: Traits class for polygon partitioning functions
// ======================================================================

#ifndef PARTITION_TRAITS_2_H
#define PARTITION_TRAITS_2_H

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_traits_2.h>
#include <CGAL/Partition_traits_2_base.h>
#include <CGAL/polygon_function_objects.h>
#include <list>

namespace CGAL {

template <class Diagonal_, class Kernel_>
class Partition_traits_2  : public Partition_traits_2_base<Kernel_>
{
  private:
    typedef Kernel_                                     Kernel;
    typedef Partition_traits_2<Diagonal_, Kernel_> Self;
  public:
    typedef Diagonal_                                   Diagonal;
    typedef CGAL::Polygon_traits_2<Kernel_>             Poly_Traits;
    typedef typename Poly_Traits::Point_2               Point_2;
    typedef ::std::list<Point_2>                        Container;
    typedef CGAL::Polygon_2<Poly_Traits, Container>     Polygon_2;
    typedef typename Kernel::Less_yx_2                  Less_yx_2;
    typedef typename Kernel::Less_xy_2                  Less_xy_2;
    typedef typename Kernel::Left_turn_2                Left_turn_2;
    typedef typename Kernel::Orientation_2              Orientation_2;
    typedef typename Kernel::Compare_y_2                Compare_y_2;
       
    typedef CGAL::Is_convex_2<Self>                     Is_convex_2;
    typedef CGAL::Is_y_monotone_2<Self>                 Is_y_monotone_2;

    // needed by Indirect_edge_compare, used in y_monotone and greene_approx
    typedef typename Kernel::Line_2                     Line_2;
    typedef typename Kernel::Construct_line_2           Construct_line_2;
    typedef typename Kernel::Compare_x_at_y_2           Compare_x_at_y_2;
    typedef typename Kernel::Is_horizontal_2            Is_horizontal_2;

    // needed by visibility graph and thus by optimal convex
    typedef typename Kernel::Ray_2                      Ray_2; 
    typedef typename Kernel::Collinear_are_ordered_along_line_2
                                            Collinear_are_ordered_along_line_2;
    typedef typename Kernel::Are_strictly_ordered_along_line_2
                                            Are_strictly_ordered_along_line_2;
    typedef typename Kernel::Intersect_2                Intersect_2;
    typedef typename Kernel::Assign_2                   Assign_2;
    typedef typename Kernel::Object_2                   Object_2;

    // needed by approx_convex (for constrained triangulation)
    // and optimal convex (for vis. graph)
    typedef typename Kernel::Segment_2                  Segment_2;
    // needed by optimal convex (for vis. graph)
    typedef typename Kernel::Construct_segment_2        Construct_segment_2;
    typedef typename Kernel::Construct_ray_2            Construct_ray_2;

 
    Construct_line_2
    construct_line_2_object() const
    {  return Construct_line_2(); }

    Compare_x_at_y_2
    compare_x_at_y_2_object() const
    { return Compare_x_at_y_2(); }

    Construct_segment_2
    construct_segment_2_object() const
    { return Construct_segment_2(); }

    Construct_ray_2
    construct_ray_2_object() const
    { return Construct_ray_2(); }

    Collinear_are_ordered_along_line_2
    collinear_are_ordered_along_line_2_object() const
    { return Collinear_are_ordered_along_line_2(); }

    Are_strictly_ordered_along_line_2
    are_strictly_ordered_along_line_2_object() const
    { return Are_strictly_ordered_along_line_2(); }

    Is_horizontal_2
    is_horizontal_2_object() const
    {  return Is_horizontal_2(); }

    Is_convex_2
    is_convex_2_object(const Self& traits) const
    {  return Is_convex_2(traits); }

    Is_y_monotone_2
    is_y_monotone_2_object(const Self& traits) const
    {  return Is_y_monotone_2(traits); }

    Intersect_2
    intersect_2_object() const
    {  return Intersect_2(); }

    Assign_2
    assign_2_object() const
    {  return Assign_2(); }

    Left_turn_2
    leftturn_2_object() const 
    { return Left_turn_2(); }

};

}

#endif // PARTITION_TRAITS_2_H
