// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.1-I-28 $
// release_date  : $CGAL_Date: 1999/10/12 $
//
// file          : include/CGAL/Cartesian.h
// package       : Cartesian_basic (3.3.6)
// source        : include/CGAL/Cartesian.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis
//
// ======================================================================


#ifndef CGAL_CARTESIAN_H
#define CGAL_CARTESIAN_H

#include <CGAL/basic.h>
#include <CGAL/cartesian_classes.h>
#include <CGAL/Cartesian_2.h>
#include <CGAL/Cartesian_3.h>
#include <CGAL/Cartesian_dynamic_d.h>

#define CGAL_REP_CLASS_DEFINED
#define CGAL_CARTESIAN_CLASS_DEFINED

CGAL_BEGIN_NAMESPACE

template< class R, class _FT >
struct Cartesian_base :
    public Cartesian_base_2<R,_FT>
    , public Cartesian_base_3<R,_FT>
    , public Cartesian_base_dynamic_d<R,_FT>
{
    // Number types and representation tag (to avoid ambiguity in
    // inheritance tree)
    typedef _FT                                           RT;
    typedef _FT                                           FT;
    typedef Cartesian_tag                                 Rep_tag;

    // All the classes are inherited, but because we inherit from a
    // template parameter, we need to explicitly write the inheritance
    // (see mail from Michael Hoffmann of July 28th 1999 in cgal-develop)

    typedef Cartesian_base_2<R,_FT>                       Kernel_base_2;
    typedef Cartesian_base_3<R,_FT>                       Kernel_base_3;
    typedef Cartesian_base_dynamic_d<R,_FT>               Kernel_base_d;

    typedef typename Kernel_base_2::Point_2               Point_2;
    typedef typename Kernel_base_2::Vector_2              Vector_2;
    typedef typename Kernel_base_2::Direction_2           Direction_2;
    typedef typename Kernel_base_2::Segment_2             Segment_2;
    typedef typename Kernel_base_2::Line_2                Line_2;
    typedef typename Kernel_base_2::Ray_2                 Ray_2;
    typedef typename Kernel_base_2::Triangle_2            Triangle_2;
    typedef typename Kernel_base_2::Circle_2              Circle_2;
    typedef typename Kernel_base_2::Iso_rectangle_2       Iso_rectangle_2;
    typedef typename Kernel_base_2::Aff_transformation_2  Aff_transformation_2;

    typedef typename Kernel_base_2::Data_accessor_2       Data_accessor_2;
    typedef typename Kernel_base_2::Conic_2               Conic_2;

    typedef typename Kernel_base_3::Point_3               Point_3;
    typedef typename Kernel_base_3::Vector_3              Vector_3;
    typedef typename Kernel_base_3::Direction_3           Direction_3;
    typedef typename Kernel_base_3::Line_3                Line_3;
    typedef typename Kernel_base_3::Plane_3               Plane_3;
    typedef typename Kernel_base_3::Ray_3                 Ray_3;
    typedef typename Kernel_base_3::Segment_3             Segment_3;
    typedef typename Kernel_base_3::Triangle_3            Triangle_3;
    typedef typename Kernel_base_3::Tetrahedron_3         Tetrahedron_3;
    typedef typename Kernel_base_3::Aff_transformation_3  Aff_transformation_3;

    typedef typename Kernel_base_d::Point_d                     Point_d;
};

CGAL_END_NAMESPACE

#include <CGAL/Kernel/Construction_objects_2.h>
#include <CGAL/Kernel/Predicate_objects_2.h>
#include <CGAL/Kernel/function_objects.h>

CGAL_BEGIN_NAMESPACE

template< class _FT >
struct Cartesian : public Cartesian_base< Cartesian<_FT>, _FT >
{
    // Number types and representation tag (to avoid ambiguity)
    typedef _FT                                           RT;
    typedef _FT                                           FT;
    typedef Cartesian_tag                                 Rep_tag;

    typedef Cartesian<FT>                                 Self;
    typedef Cartesian_base<Self,FT>                       Kernel_base;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
    // The other classes are inherited and because of partial specialization,
    // Cartesian<FT>::Point_2 is exactly CGAL::Point_2< Cartesian<FT> >
    // As above, we still need to write down the inheritance explicitly

    typedef typename Kernel_base::Point_2                 Point_2;
    typedef typename Kernel_base::Vector_2                Vector_2;
    typedef typename Kernel_base::Direction_2             Direction_2;
    typedef typename Kernel_base::Segment_2               Segment_2;
    typedef typename Kernel_base::Line_2                  Line_2;
    typedef typename Kernel_base::Ray_2                   Ray_2;
    typedef typename Kernel_base::Triangle_2              Triangle_2;
    typedef typename Kernel_base::Circle_2                Circle_2;
    typedef typename Kernel_base::Iso_rectangle_2         Iso_rectangle_2;
    typedef typename Kernel_base::Aff_transformation_2    Aff_transformation_2;

    typedef typename Kernel_base::Data_accessor_2         Data_accessor_2;
    typedef typename Kernel_base::Conic_2                 Conic_2;

    typedef typename Kernel_base::Point_3                 Point_3;
    typedef typename Kernel_base::Vector_3                Vector_3;
    typedef typename Kernel_base::Direction_3             Direction_3;
    typedef typename Kernel_base::Line_3                  Line_3;
    typedef typename Kernel_base::Plane_3                 Plane_3;
    typedef typename Kernel_base::Ray_3                   Ray_3;
    typedef typename Kernel_base::Segment_3               Segment_3;
    typedef typename Kernel_base::Triangle_3              Triangle_3;
    typedef typename Kernel_base::Tetrahedron_3           Tetrahedron_3;
    typedef typename Kernel_base::Aff_transformation_3    Aff_transformation_3;

    typedef typename Kernel_base::Point_d                 Point_d;
    typedef typename Kernel_base::Vector_d                Vector_d;
    typedef typename Kernel_base::Direction_d             Direction_d;

#else
    // Now CGAL::Point_2<R> is only a wrapper around CGAL::PointC2<R>
    // It is necessary to redefine here the classes to ensure that
    // Cartesian<FT>::Point_2 is exactly CGAL::Point_2< Cartesian<FT> >

    // Cartesian<FT>::Base is needed so that CGAL::Point_2< Cartesian<FT> >
    // can inherit from Cartesian<FT>::Point_2_base

    typedef typename Kernel_base::Point_2                 Point_2_base;
    typedef typename Kernel_base::Vector_2                Vector_2_base;
    typedef typename Kernel_base::Direction_2             Direction_2_base;
    typedef typename Kernel_base::Segment_2               Segment_2_base;
    typedef typename Kernel_base::Line_2                  Line_2_base;
    typedef typename Kernel_base::Ray_2                   Ray_2_base;
    typedef typename Kernel_base::Triangle_2              Triangle_2_base;
    typedef typename Kernel_base::Circle_2                Circle_2_base;
    typedef typename Kernel_base::Iso_rectangle_2         Iso_rectangle_2_base;
    typedef typename Kernel_base::Aff_transformation_2    Aff_transformation_2_base;

    typedef typename Kernel_base::Point_3                 Point_3_base;
    typedef typename Kernel_base::Vector_3                Vector_3_base;
    typedef typename Kernel_base::Direction_3             Direction_3_base;
    typedef typename Kernel_base::Line_3                  Line_3_base;
    typedef typename Kernel_base::Plane_3                 Plane_3_base;
    typedef typename Kernel_base::Ray_3                   Ray_3_base;
    typedef typename Kernel_base::Segment_3               Segment_3_base;
    typedef typename Kernel_base::Triangle_3              Triangle_3_base;
    typedef typename Kernel_base::Tetrahedron_3           Tetrahedron_3_base;
    typedef typename Kernel_base::Aff_transformation_3    Aff_transformation_3_base;
  
    typedef typename Kernel_base::Point_d                 Point_d_base;
    typedef typename Kernel_base::Vector_d                Vector_d_base;
    typedef typename Kernel_base::Direction_d             Direction_d_base;

    // Note: necessary to qualify Point_2 by CGAL:: to disambiguate between
    // Point_2 in the current namespace (nested within CGAL) and
    // CGAL::Point_2< Cartesian<FT> > (which is in the CGAL namespace)

    typedef CGAL::Point_2<Self>                           Point_2;
    typedef CGAL::Vector_2<Self>                          Vector_2;
    typedef CGAL::Direction_2<Self>                       Direction_2;
    typedef CGAL::Line_2<Self>                            Line_2;
    typedef CGAL::Ray_2<Self>                             Ray_2;
    typedef CGAL::Segment_2<Self>                         Segment_2;
    typedef CGAL::Triangle_2<Self>                        Triangle_2;
    typedef CGAL::Circle_2<Self>                          Circle_2;
    typedef CGAL::Iso_rectangle_2<Self>                   Iso_rectangle_2;
    typedef CGAL::Aff_transformation_2<Self>              Aff_transformation_2;

    typedef Data_accessorC2<Self>                         Data_accessor_2;
    typedef ConicCPA2<Point_2,Data_accessor_2>            Conic_2;

    typedef CGAL::Point_3<Self>                           Point_3;
    typedef CGAL::Vector_3<Self>                          Vector_3;
    typedef CGAL::Direction_3<Self>                       Direction_3;
    typedef CGAL::Line_3<Self>                            Line_3;
    typedef CGAL::Plane_3<Self>                           Plane_3;
    typedef CGAL::Ray_3<Self>                             Ray_3;
    typedef CGAL::Segment_3<Self>                         Segment_3;
    typedef CGAL::Triangle_3<Self>                        Triangle_3;
    typedef CGAL::Tetrahedron_3<Self>                     Tetrahedron_3;
    typedef CGAL::Aff_transformation_3<Self>              Aff_transformation_3;

    typedef CGAL::Point_d<Self>                           Point_d;
    typedef CGAL::Vector_d<Self>                          Vector_d;
    typedef CGAL::Director_d<Self>                        Director_d;

#endif // CGAL_CFG_NO_ADVANCED_KERNEL

    // TODO: cleanup
    static   FT make_FT(const RT & num, const RT& denom) { return num/denom;}
    static   FT make_FT(const RT & num)                  { return num;}
    static   RT FT_numerator(const FT &r)                { return r;}
    static   RT FT_denominator(const FT &)               { return RT(1);}

typedef CGALi::Construct<Point_2>              Construct_point_2;
typedef CGALi::Construct<Vector_2>             Construct_vector_2;
typedef CGALi::Construct<Direction_2>          Construct_direction_2;
typedef CGALi::Construct<Segment_2>            Construct_segment_2;
typedef CGALi::Construct<Line_2>               Construct_line_2;
typedef CGALi::Construct<Ray_2>                Construct_ray_2;
typedef CGALi::Construct<Circle_2>             Construct_circle_2;
typedef CGALi::Construct<Triangle_2>           Construct_triangle_2;
typedef CGALi::Construct<Aff_transformation_2> Construct_aff_transformation_2;

Construct_point_2 
construct_point_2_object() const 
{ return Construct_point_2(); }

Construct_vector_2
construct_vector_2_object() const 
{ return Construct_vector_2(); }

Construct_direction_2
construct_direction_2_object() const 
{ return Construct_direction_2(); }

Construct_segment_2
construct_segment_2_object() const 
{ return Construct_segment_2(); }

Construct_line_2
construct_line_2_object() const 
{ return Construct_line_2(); }

Construct_ray_2
construct_ray_2_object() const 
{ return Construct_ray_2(); }

Construct_circle_2
construct_circle_2_object() const 
{ return Construct_circle_2(); }

Construct_triangle_2
construct_triangle_2_object() const 
{ return Construct_triangle_2(); }

Construct_aff_transformation_2
construct_aff_transformation_2_object() const 
{ return Construct_aff_transformation_2(); }

typedef CGALi::Call_point_to_get<Point_2>              Construct_point_on_2;
Construct_point_on_2
construct_point_on_2_object() const 
{ return Construct_point_on_2(); }

typedef CGALi::Call_second_point_to_get<Point_2>  Construct_second_point_on_2;
Construct_second_point_on_2
construct_second_point_on_2_object() const 
{ return Construct_second_point_on_2(); }

typedef CGALi::Call_source_to_get<Point_2>        Construct_source_point_2;
Construct_source_point_2
construct_source_point_2_object() const 
{ return Construct_source_point_2(); }

typedef CGALi::Call_target_to_get<Point_2>        Construct_target_point_2;
Construct_target_point_2
construct_target_point_2_object() const 
{ return Construct_target_point_2(); }

typedef CGALi::Call_min_to_get<Point_2>           Construct_min_point_2;
Construct_min_point_2
construct_min_point_2_object() const 
{ return Construct_min_point_2(); }

typedef CGALi::Call_max_to_get<Point_2>           Construct_max_point_2;
Construct_max_point_2
construct_max_point_2_object() const 
{ return Construct_max_point_2(); }

typedef CGALi::Call_direction_to_get<Direction_2> Construct_direction_of_line_2;
Construct_direction_of_line_2
construct_direction_of_line_2_object() const 
{ return Construct_direction_of_line_2(); }

typedef CGALi::Call_direction_to_get<Direction_2> Construct_direction_of_ray_2;
Construct_direction_of_ray_2
construct_direction_of_ray_2_object() const 
{ return Construct_direction_of_ray_2(); }

typedef CGALi::Call_supporting_line_to_get<Line_2> Construct_supporting_line_2;
Construct_supporting_line_2
construct_supporting_line_2_object() const 
{ return Construct_supporting_line_2(); }

typedef CGALi::Call_perpendicular_to_get<Vector_2> Construct_perpendicular_vector_2;
Construct_perpendicular_vector_2
construct_perpendicular_vector_2_object() const 
{ return Construct_perpendicular_vector_2(); }

typedef CGALi::Call_perpendicular_to_get<Direction_2>  Construct_perpendicular_direction_2;
Construct_perpendicular_direction_2
construct_perpendicular_direction_2_object() const 
{ return Construct_perpendicular_direction_2(); }

typedef CGALi::Call_perpendicular_to_get<Line_2>   Construct_perpendicular_line_2;
Construct_perpendicular_line_2
construct_perpendicular_line_2_object() const 
{ return Construct_perpendicular_line_2(); }

typedef CGALi::p_Midpoint<Point_2>                 Construct_midpoint;
Construct_midpoint
construct_midpoint_object() const 
{ return Construct_midpoint(); }

typedef CGALi::p_Circumcenter<Point_2>             Construct_circumcenter_2;
Construct_circumcenter_2
construct_circumcenter_2_object() const 
{ return Construct_circumcenter_2(); }

typedef CGALi::pl_Bisector<Point_2, Line_2>        Construct_bisector_2;
Construct_bisector_2
construct_bisector_2_object() const 
{ return Construct_bisector_2(); }

typedef CGALi::Call_opposite_to_get<Segment_2>     Construct_opposite_segment_2;
Construct_opposite_segment_2
construct_opposite_segment_2_object() const 
{ return Construct_opposite_segment_2(); }

typedef CGALi::Call_opposite_to_get<Ray_2>         Construct_opposite_ray_2;
Construct_opposite_ray_2
construct_opposite_ray_2_object() const 
{ return Construct_opposite_ray_2(); }

typedef CGALi::Call_opposite_to_get<Line_2>        Construct_opposite_line_2;
Construct_opposite_line_2
construct_opposite_line_2_object() const 
{ return Construct_opposite_line_2(); }

typedef CGALi::Call_opposite_to_get<Triangle_2>    Construct_opposite_triangle_2;
Construct_opposite_triangle_2
construct_opposite_triangle_2_object() const 
{ return Construct_opposite_triangle_2(); }

typedef CGALi::Call_opposite_to_get<Circle_2>      Construct_opposite_circle_2;
Construct_opposite_circle_2
construct_opposite_circle_2_object() const 
{ return Construct_opposite_circle_2(); }

typedef CGALi::Call_transform                      Transform_2;
Transform_2
transform_2_object() const 
{ return Transform_2(); }

typedef CGALi::Intersect                           Intersect_2;
Intersect_2
intersect_2_object() const 
{ return Intersect_2(); }

typedef CGALi::Call_y_at_x_to_get<FT>              Compute_y_at_x;
Compute_y_at_x
compute_y_at_x_object() const 
{ return Compute_y_at_x(); }

typedef CGALi::Call_squared_length_to_get<FT>      Compute_squared_length_2;
Compute_squared_length_2
Compute_squared_length_2_object() const 
{ return Compute_squared_length_2(); }

typedef CGALi::Equal                               Equal_2;
Equal_2
equal_2_object() const 
{ return Equal_2(); }

typedef CGALi::Equal_x                             Equal_x_2;
Equal_x_2
equal_x_2_object() const 
{ return Equal_x_2(); }

typedef CGALi::Equal_y                             Equal_y_2;
Equal_y_2
equal_y_2_object() const 
{ return Equal_y_2(); }

typedef CGALi::Less_x                              Less_x_2;
Less_x_2
less_x_2_object() const 
{ return Less_x_2(); }

typedef CGALi::Less_y                              Less_y_2;
Less_y_2
less_y_2_object() const 
{ return Less_y_2(); }

typedef CGAL::p_Less_xy<Point_2>                   Less_xy_2;
Less_xy_2
less_xy_2_object() const 
{ return Less_xy_2(); }

typedef CGALi::Compare_x                           Compare_x_2;
Compare_x_2
compare_x_2_object() const 
{ return Compare_x_2(); }

typedef CGALi::Compare_y                           Compare_y_2;
Compare_y_2
compare_y_2_object() const 
{ return Compare_y_2(); }

typedef CGALi::Compare_xy                          Compare_xy_2;
Compare_xy_2
compare_xy_2_object() const 
{ return Compare_xy_2(); }

typedef CGALi::Compare_y_at_x                      Compare_y_at_x_2;
Compare_y_at_x_2
compare_y_at_x_2_object() const 
{ return Compare_y_at_x_2(); }

typedef CGAL ::p_Less_dist_to_line_2p<Point_2>     Less_signed_distance_to_line_2;
Less_signed_distance_to_line_2
less_signed_distance_to_line_2_object(const Point_2& p, const Point_2& q) const 
{ return Less_signed_distance_to_line_2(p,q); }

typedef CGAL ::p_Less_rotate_ccw<Point_2>          Less_rotate_ccw_2;
Less_rotate_ccw_2
less_rotate_ccw_2(const Point_2& p) const 
{ return Less_rotate_ccw_2(p); }

typedef CGAL ::p_Leftturn<Point_2>                 Leftturn_2;
Leftturn_2
leftturn_2_object() const 
{ return Leftturn_2(); }

typedef CGAL ::p_Left_of_line_2p<Point_2>          Left_of_line_2;
Left_of_line_2
left_of_line_2_object(const Point_2& p, const Point_2& q) const 
{ return Left_of_line_2(p,q); }

typedef CGALi::Collinear                           Collinear_2;
Collinear_2
collinear_2_object() const 
{ return Collinear_2(); }

typedef CGAL ::p_Orientation<Point_2>              Orientation_2;
Orientation_2
orientation_2_object() const 
{ return Orientation_2(); }

typedef CGALi::Side_of_oriented_circle             Side_of_oriented_circle_2;
Side_of_oriented_circle_2
side_of_oriented_circle_2_object() const 
{ return Side_of_oriented_circle_2(); }

typedef CGALi::Side_of_bounded_circle              Side_of_bounded_circle_2;
Side_of_bounded_circle_2
side_of_bounded_circle_2_object() const
{ return Side_of_bounded_circle_2(); }

typedef CGALi::Call_is_horizontal                  Is_horizontal_2;
Is_horizontal_2
is_horizontal_2_object() const 
{ return Is_horizontal_2(); }

typedef CGALi::Call_is_vertical                    Is_vertical_2;
Is_vertical_2
is_vertical_2_object() const 
{ return Is_vertical_2(); }

typedef CGALi::Call_is_degenerate                  Is_degenerate_2;
Is_degenerate_2
is_degenerate_2_object() const 
{ return Is_degenerate_2(); }

typedef CGALi::Call_has_on                         Has_on_2;
Has_on_2
has_on_2_object() const 
{ return Has_on_2(); }

typedef CGALi::Call_collinear_has_on               Collinear_has_on_2;
Collinear_has_on_2
collinear_has_on_2_object() const 
{ return Collinear_has_on_2(); }

typedef CGALi::Call_has_on_bounded_side            Has_on_bounded_side_2;
Has_on_bounded_side_2
has_on_bounded_side_2_object() const 
{ return Has_on_bounded_side_2(); }

typedef CGALi::Call_has_on_unbounded_side          Has_on_unbounded_side_2;
Has_on_unbounded_side_2
has_on_unbounded_side_2_object() const 
{ return Has_on_unbounded_side_2(); }

typedef CGALi::Call_has_on_boundary                Has_on_boundary_2;
Has_on_boundary_2
has_on_boundary_2_object() const 
{ return Has_on_boundary_2(); }

typedef CGALi::Call_has_on_positive_side           Has_on_positive_side_2;
Has_on_positive_side_2
has_on_positive_side_2_object() const 
{ return Has_on_positive_side_2(); }

typedef CGALi::Call_has_on_negative_side           Has_on_negative_side_2;
Has_on_negative_side_2
has_on_negative_side_2_object() const 
{ return Has_on_negative_side_2(); }

typedef CGALi::Call_oriented_side                  Oriented_side_2;
Oriented_side_2
oriented_side_2_object() const 
{ return Oriented_side_2(); }

typedef CGALi::Are_ordered_along_line              Are_ordered_along_line_2 ;
Are_ordered_along_line_2
are_ordered_along_line_2_object() const 
{ return Are_ordered_along_line_2(); }

typedef CGALi::Are_strictly_ordered_along_line     Are_strictly_ordered_along_line_2;
Are_strictly_ordered_along_line_2
are_strictly_ordered_along_line_2_object() const 
{ return Are_strictly_ordered_along_line_2(); }

typedef CGALi::Collinear_are_ordered_along_line    Collinear_are_ordered_along_line_2;
Collinear_are_ordered_along_line_2
collinear_are_ordered_along_line_2_object() const 
{ return Collinear_are_ordered_along_line_2(); }

typedef CGALi::Collinear_are_strictly_ordered_along_line Collinear_are_strictly_ordered_along_line_2;
Collinear_are_strictly_ordered_along_line_2
collinear_are_strictly_ordered_along_line_2_object() const 
{ return Collinear_are_strictly_ordered_along_line_2(); }

typedef CGALi::Construct<Point_3>                  Construct_point_3;
typedef CGALi::Construct<Vector_3>                 Construct_vector_3;
typedef CGALi::Construct<Direction_3>              Construct_direction_3;
typedef CGALi::Construct<Segment_3>                Construct_segment_3;
typedef CGALi::Construct<Plane_3>                  Construct_plane_3;
typedef CGALi::Construct<Line_3>                   Construct_line_3;
typedef CGALi::Construct<Ray_3>                    Construct_ray_3;
typedef CGALi::Construct<Triangle_3>               Construct_triangle_3;
typedef CGALi::Construct<Tetrahedron_3>            Construct_tetrahedron_3;
typedef CGALi::Construct<Aff_transformation_3>     Construct_aff_transformation_3;

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
construct_tetrahedron_3_object() const 
{ return Construct_tetrahedron_3(); }

Construct_aff_transformation_3
construct_aff_transformation_3_object() const 
{ return Construct_aff_transformation_3(); }

typedef CGALi::Call_point_to_get<Point_3>              Construct_point_on_3;
Construct_point_on_3
construct_point_on_3_object() const 
{ return Construct_point_on_3(); }

typedef CGALi::Call_second_point_to_get<Point_3>       Construct_second_point_on_3;
Construct_second_point_on_3
construct_second_point_on_3_object() const 
{ return Construct_second_point_on_3(); }

typedef CGALi::Call_perpendicular_plane_to_get<Plane_3> Construct_perpendicular_plane_3;
Construct_perpendicular_plane_3
construct_perpendicular_plane_3() const 
{ return Construct_perpendicular_plane_3(); }

typedef CGALi::p_Midpoint<Point_3>                     Construct_midpoint_3;
Construct_midpoint_3
construct_midpoint_3_object() const 
{ return Construct_midpoint_3(); }

typedef CGALi::Call_opposite_to_get<Segment_3>         Construct_opposite_segment_3;
Construct_opposite_segment_3
construct_opposite_segment_3_object() const 
{ return Construct_opposite_segment_3(); }

typedef CGALi::Call_opposite_to_get<Ray_3>             Construct_opposite_ray_3;
Construct_opposite_ray_3
construct_opposite_ray_3_object() const 
{ return Construct_opposite_ray_3(); }

typedef CGALi::Call_opposite_to_get<Line_3>            Construct_opposite_line_3;
Construct_opposite_line_3
construct_opposite_line_3_object() const 
{ return Construct_opposite_line_3(); }

typedef CGALi::Call_supporting_plane_to_get<Plane_3>   Construct_supporting_plane_3;
Construct_supporting_plane_3
construct_supporting_plane_3_object() const 
{ return Construct_supporting_plane_3(); }

typedef CGALi::Call_transform                          Transform_3;
Transform_3
transform_3_object() const 
{ return Transform_2(); }

typedef CGALi::Intersect                               Intersect_3;
Intersect_3
intersect_3_object() const 
{ return Intersect_3(); }

typedef CGALi::Call_squared_length_to_get<FT>          Compute_squared_length_3;
Compute_squared_length_3
compute_squared_length_3_object() const 
{ return Compute_squared_length_3(); }

typedef CGALi::Equal                                   Equal_3;
Equal_3
equal_3_object() const 
{ return Equal_3(); }

typedef CGALi::Equal_x                                 Equal_x_3;
Equal_x_3
equal_x_3_object() const 
{ return Equal_x_3(); }

typedef CGALi::Equal_y                                 Equal_y_3;
Equal_y_3
equal_y_3_object() const 
{ return Equal_y_3(); }

typedef CGALi::Equal_z                                 Equal_z_3;
Equal_z_3
equal_z_3_object() const 
{ return Equal_z_3(); }

typedef CGALi::Equal_xy                                Equal_xy_3;
Equal_xy_3
equal_xy_3_object() const 
{ return Equal_xy_3(); }

typedef CGALi::Less_x                                  Less_x_3;
Less_x_3
less_x_3_object() const 
{ return Less_x_3(); }

typedef CGALi::Less_y                                  Less_y_3;
Less_y_3
less_y_3_object() const 
{ return Less_y_3(); }

typedef CGALi::Less_z                                  Less_z_3;
Less_z_3
less_z_3_object() const 
{ return Less_z_3(); }

typedef CGAL::p_Less_xy<Point_3>                       Less_xy_3;
Less_xy_3
less_xy_3_object() const 
{ return Less_xy_3(); }

typedef CGALi::Less_xyz                                Less_xyz_3;
Less_xyz_3
less_xyz_3_object() const 
{ return Less_xyz_3(); }

typedef CGALi::Compare_x                               Compare_x_3;
Compare_x_3
compare_x_3_object() const 
{ return Compare_x_3(); }

typedef CGALi::Compare_y                               Compare_y_3;
Compare_y_3
compare_y_3_object() const 
{ return Compare_y_3(); }

typedef CGALi::Compare_z                               Compare_z_3;
Compare_z_3
compare_z_3_object() const 
{ return Compare_z_3(); }

typedef CGALi::Compare_xy                              Compare_xy_3;
Compare_xy_3
compare_xy_3_object() const 
{ return Compare_xyz_3(); }

typedef CGALi::Compare_xyz                             Compare_xyz_3;
Compare_xyz_3
compare_xyz_3_object() const 
{ return Compare_xyz_3(); }

typedef CGALi::Collinear                               Collinear_3;
Collinear_3
collinear_3_object() const 
{ return Collinear_3(); }

typedef CGALi::Coplanar                                Coplanar_3 ;
Coplanar_3
coplanar_3_object() const 
{ return Coplanar_3(); }

typedef CGAL ::p_Orientation<Point_3>                  Orientation_3;
Orientation_3
orientation_3_object() const 
{ return Orientation_3(); }

typedef CGALi::Call_is_degenerate                      Is_degenerate_3;
Is_degenerate_3
is_degenerate_3_object() const 
{ return Is_degenerate_3(); }

typedef CGALi::Call_has_on                             Has_on_3;
Has_on_3
has_on_3_object() const 
{ return Has_on_3(); }

typedef CGALi::Call_has_on_bounded_side                Has_on_bounded_side_3;
Has_on_bounded_side_3
has_on_bounded_side_3_object() const 
{ return Has_on_bounded_side_3(); }

typedef CGALi::Call_has_on_unbounded_side              Has_on_unbounded_side_3;
Has_on_unbounded_side_3
has_on_unbounded_side_3_object() const 
{ return Has_on_unbounded_side_3(); }

typedef CGALi::Call_has_on_boundary                    Has_on_boundary_3;
Has_on_boundary_3
has_on_boundary_3_object() const 
{ return Has_on_boundary_3(); }

typedef CGALi::Call_has_on_positive_side               Has_on_positive_side_3;
Has_on_positive_side_3
has_on_positive_side_3_object() const 
{ return Has_on_positive_side_3(); }

typedef CGALi::Call_has_on_negative_side               Has_on_negative_side_3;
Has_on_negative_side_3
has_on_negative_side_3_object() const 
{ return Has_on_negative_side_3(); }

typedef CGALi::Call_oriented_side                      Oriented_side_3;
Oriented_side_3
oriented_side_3_object() const 
{ return Oriented_side_3(); }

typedef CGALi::Are_ordered_along_line                  Are_ordered_along_line_3 ;
Are_ordered_along_line_3
are_ordered_along_line_3_object() const 
{ return Are_ordered_along_line_3(); }

typedef CGALi::Are_strictly_ordered_along_line         Are_strictly_ordered_along_line_3;
Are_strictly_ordered_along_line_3
are_strictly_ordered_along_line_3_object() const 
{ return Are_strictly_ordered_along_line_3(); }

typedef CGALi::Collinear_are_ordered_along_line        Collinear_are_ordered_along_line_3;
Collinear_are_ordered_along_line_3
collinear_are_ordered_along_line_3_object() const 
{ return Collinear_are_ordered_along_line_3(); }

typedef CGALi::Collinear_are_strictly_ordered_along_line Collinear_are_strictly_ordered_along_line_3;
Collinear_are_strictly_ordered_along_line_3
collinear_are_strictly_ordered_along_line_3_object() const 
{ return Collinear_are_strictly_ordered_along_line_3(); }

typedef CGALi::Side_of_oriented_sphere                 Side_of_oriented_sphere_3;
Side_of_oriented_sphere_3
side_of_oriented_sphere_3_object() const 
{ return Side_of_oriented_sphere_3(); }

typedef CGALi::Side_of_bounded_sphere                  Side_of_bounded_sphere_3;
Side_of_bounded_sphere_3
side_of_bounded_sphere_3_object() const 
{ return Side_of_bounded_sphere_3(); }

};

CGAL_END_NAMESPACE

#include <CGAL/iterator_traits_pointer_specs_for_cartesian_kernel.h>

#endif // CGAL_CARTESIAN_H
