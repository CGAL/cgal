// ============================================================================
//
// Copyright (c) 1997  The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Triangulation_euclidean_traits_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
#define CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/basic_constructions_2.h>
#include <CGAL/distance_predicates_2.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>

CGAL_BEGIN_NAMESPACE 

template < class R >
class Triangulation_euclidean_traits_2 {
public:
  typedef R Rep;
  typedef Point_2<R>  Point_2;
  typedef Segment_2<R> Segment_2;
  typedef Triangle_2<R> Triangle_2;
  typedef Line_2<R> Line_2;
  typedef Direction_2<R> Direction_2;
  typedef Ray_2<R> Ray_2;

  typedef typename R::Compare_x_2                Compare_x_2;
  typedef typename R::Compare_y_2                Compare_y_2;
  typedef typename R::Orientation_2              Orientation_2;
  typedef typename R::Side_of_oriented_circle_2  Side_of_oriented_circle_2;
  typedef typename R::Construct_circumcenter_2   Construct_circumcenter_2;
  typedef typename R::Construct_bisector_2       Construct_bisector_2;
  //typedef typename R::Construct_midpoint         Construct_midpoint;
  typedef typename R::Compare_distance_2         Compare_distance_2;
  typedef typename R::Construct_segment_2        Construct_segment_2;
  typedef typename R::Construct_triangle_2       Construct_triangle_2;
  //typedef typename R::Construct_direction_2      Construct_direction_2;
  typedef typename R::Construct_ray_2            Construct_ray_2;
  typedef typename R::Construct_direction_of_line_2
                                           Construct_direction_of_line_2;

  // for compatibility with previous versions
  typedef Point_2      Point;
  typedef Segment_2    Segment;
  typedef Triangle_2   Triangle;
  typedef Ray_2        Ray;
  typedef Line_2       Line;
  typedef Direction_2  Direction;

  Triangulation_euclidean_traits_2() {}
  Triangulation_euclidean_traits_2(const Triangulation_euclidean_traits_2 &) {}
  Triangulation_euclidean_traits_2 &operator=
      (const Triangulation_euclidean_traits_2 &)
  {return *this;}
 
  Compare_x_2
  compare_x_2_object() const
    { return Compare_x_2();}

  Compare_y_2
  compare_y_2_object() const
    { return Compare_y_2();}
  
  Orientation_2
  orientation_2_object() const
    { return Orientation_2();}

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
    {return Side_of_oriented_circle_2();}
 
  Construct_circumcenter_2
  construct_circumcenter_2_object() const
    { return Construct_circumcenter_2();}

  Construct_bisector_2
  construct_bisector_2_object() const
    {return Construct_bisector_2();}
  
 //  Construct_midpoint
//   construct_midpoint_object() const
//     {return Construct_midpoint();}


  Compare_distance_2
  compare_distance_2_object() const
    {return Compare_distance_2();}

  Construct_direction_of_line_2
  construct_direction_of_line_2_object() const
    {return  Construct_direction_of_line_2();}

  Construct_segment_2  construct_segment_2_object() const
    {return Construct_segment_2();}

  Construct_triangle_2  construct_triangle_2_object() const
    {return Construct_triangle_2();}

//   Construct_direction_2  construct_direction_2_object() const
//     {return Construct_direction_2();}

  Construct_ray_2  construct_ray_2_object() const
    {return Construct_ray_2();}

};

CGAL_END_NAMESPACE 

#endif // CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
