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
// release       :
// release_date  :
//
// file          : include/CGAL/Kernel/Construction_objects_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_KERNEL_CONSTRUCTION_OBJECTS_3_H
#define CGAL_KERNEL_CONSTRUCTION_OBJECTS_3_H

#include <CGAL/Kernel/function_objects.h>

CGAL_BEGIN_NAMESPACE

template < class R >
class Kernel_construction_objects_3
{
public:
typedef typename R::FT                         FT;
typedef typename R::RT                         RT;
typedef typename R::Point_3                    Point_3;
typedef typename R::Vector_3                   Vector_3;
typedef typename R::Direction_3                Direction_3;
typedef typename R::Line_3                     Line_3;
typedef typename R::Plane_3                    Plane_3;
typedef typename R::Ray_3                      Ray_3;
typedef typename R::Segment_3                  Segment_3;
typedef typename R::Triangle_3                 Triangle_3;
typedef typename R::Tetrahedron_3              Tetrahedron_3;
typedef typename R::Aff_transformation_3       Aff_transformation_3;

typedef CGALi::Construct<Point_3>              Construct_point_3;
typedef CGALi::Construct<Vector_3>             Construct_vector_3;
typedef CGALi::Construct<Direction_3>          Construct_direction_3;
typedef CGALi::Construct<Segment_3>            Construct_segment_3;
typedef CGALi::Construct<Plane_3>              Construct_plane_3;
typedef CGALi::Construct<Line_3>               Construct_line_3;
typedef CGALi::Construct<Ray_3>                Construct_ray_3;
typedef CGALi::Construct<Triangle_3>           Construct_triangle_3;
typedef CGALi::Construct<Tetrahedron_3>        Construct_tetrahedron_3;
typedef CGALi::Construct<Aff_transformation_3> Construct_aff_transformation_3;

Construct_point_3 
construct_point_3_object() const 
{ return Construct_point_3(); }

Construct_vector_3
construct_vector_3_object() const 
{ return Construct_vector_3(); }

Construct_direction_3
construct_direction_3_object() const 
{ return Construct_direction_3(); }

Construct_segment_3
construct_segment_3_object() const 
{ return Construct_segment_3(); }

Construct_plane_3
construct_plane_3_object() const 
{ return Construct_plane_3(); }

Construct_line_3
construct_line_3_object() const 
{ return Construct_line_3(); }

Construct_ray_3
construct_ray_3_object() const 
{ return Construct_ray_3(); }

Construct_triangle_3
construct_triangle_3_object() const 
{ return Construct_triangle_3(); }

Construct_tetrahedron_3
construct_tetrahedron_object() const 
{ return Construct_tetrahedron_3(); }

Construct_aff_transformation_3
construct_aff_transformation_3_object() const 
{ return Construct_aff_transformation_3(); }

typedef CGALi::Call_point_to_get<Point_3>              Construct_point_on_3;
Construct_point_on_3
construct_point_on_3_object() const 
{ return Construct_point_on_3(); }

typedef CGALi::Call_second_point_to_get<Point_3>  Construct_second_point_on_3;
Construct_second_point_on_3
construct_second_point_on_3_object() const 
{ return Construct_second_point_on_3(); }

typedef CGALi::Call_perpendicular_plane_to_get<Plane_3> 
                                         Construct_perpendicular_plane_3;
Construct_perpendicular_plane_3
construct_perpendicular_plane_3() const 
{ return Construct_perpendicular_plane_3(); }

typedef CGALi::p_Midpoint<Point_3>                     Construct_midpoint_3;
Construct_midpoint_3
construct_midpoint_3_object() const 
{ return Construct_midpoint_3(); }

typedef CGALi::p_Circumcenter<Point_3>                Construct_circumcenter_3;
Construct_circumcenter_3
construct_circumcenter_3_object() const 
{ return Construct_circumcenter_3(); }

typedef CGALi::Call_opposite_to_get<Segment_3>  Construct_opposite_segment_3;
Construct_opposite_segment_3
construct_opposite_segment_3_object() const 
{ return Construct_opposite_segment_3(); }

typedef CGALi::Call_opposite_to_get<Ray_3>         Construct_opposite_ray_3;
Construct_opposite_ray_3
construct_opposite_ray_3_object() const 
{ return Construct_opposite_ray_3(); }

typedef CGALi::Call_opposite_to_get<Line_3>        Construct_opposite_line_3;
Construct_opposite_line_3
construct_opposite_line_3_object() const 
{ return Construct_opposite_line_3(); }

typedef CGALi::Call_supporting_plane_to_get<Plane_3> 
                                             Construct_supporting_plane_3;
Construct_supporting_plane_3
construct_supporting_plane_3_object() const 
{ return Construct_supporting_plane_3(); }

typedef CGALi::Call_transform                      Transform_3;
Transform_3
transform_3_object() const 
{ return Transform_2(); }

typedef CGALi::Intersect                           Intersect_3;
Intersect_3
intersect_3_object() const 
{ return Intersect_3(); }

typedef CGALi::Call_squared_length_to_get<FT>      Compute_squared_length_3;
Compute_squared_length_3
compute_squared_length_3_object() const 
{ return Compute_squared_length_3(); }

};

CGAL_END_NAMESPACE

// This macro is provided for convenience in defining the Kernel
// function objects inside a new representation class.
// See Cartesian_3.h and Cartesian.h

#define CGAL_UNPACK_KERNEL_CONSTRUCTION_OBJECTS_3(CO)

#endif // CGAL_KERNEL_CONSTRUCTION_OBJECTS_3_H
