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
// file          : include/CGAL/Cartesian_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_2_H
#define CGAL_CARTESIAN_2_H

#include <CGAL/basic.h>
#include <CGAL/cartesian_classes.h>

#ifdef CGAL_CFG_NO_ADVANCED_KERNEL
  // Because we cannot use Michael's scheme, we need the wrapper classes
  // We include them (they are common to Cartesian and Homogeneous)
#include <CGAL/user_classes.h>
#endif // CGAL_CFG_NO_ADVANCED_KERNEL

#define CGAL_REP_CLASS_DEFINED
#define CGAL_CARTESIAN_CLASS_DEFINED

CGAL_BEGIN_NAMESPACE

template< class R, class FT_ >
struct Cartesian_base_2
{
    // Number types and representation tag
    typedef FT_                                    RT;
    typedef FT_                                    FT;
    typedef Cartesian_tag                          Rep_tag;
    typedef Cartesian_tag                          Kernel_tag;
    typedef CGAL::Object                           Object_2;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
    typedef CGAL::Point_2<R,Rep_tag>               Point_2;
    typedef CGAL::Vector_2<R,Rep_tag>              Vector_2;
    typedef CGAL::Direction_2<R,Rep_tag>           Direction_2;
    typedef CGAL::Segment_2<R,Rep_tag>             Segment_2;
    typedef CGAL::Line_2<R,Rep_tag>                Line_2;
    typedef CGAL::Ray_2<R,Rep_tag>                 Ray_2;
    typedef CGAL::Triangle_2<R,Rep_tag>            Triangle_2;
    typedef CGAL::Circle_2<R,Rep_tag>              Circle_2;
    typedef CGAL::Iso_rectangle_2<R,Rep_tag>       Iso_rectangle_2;
    typedef CGAL::Aff_transformation_2<R,Rep_tag>  Aff_transformation_2;
    typedef CGAL::Data_accessor_2<R,Rep_tag>       Data_accessor_2;
    typedef CGAL::ConicCPA2<Point_2,Data_accessor_2> Conic_2;
#else
    typedef PointC2<R>                             Point_2;
    typedef VectorC2<R>                            Vector_2;
    typedef DirectionC2<R>                         Direction_2;
    typedef SegmentC2<R>                           Segment_2;
    typedef LineC2<R>                              Line_2;
    typedef RayC2<R>                               Ray_2;
    typedef TriangleC2<R>                          Triangle_2;
    typedef CircleC2<R>                            Circle_2;
    typedef Iso_rectangleC2<R>                     Iso_rectangle_2;
    typedef Aff_transformationC2<R>                Aff_transformation_2;
    typedef Data_accessorC2<R>                     Data_accessor_2;
    typedef ConicCPA2<Point_2,Data_accessor_2>     Conic_2;
#endif // CGAL_CFG_NO_ADVANCED_KERNEL

};

CGAL_END_NAMESPACE

#include <CGAL/Cartesian/Point_2.h>
#include <CGAL/Cartesian/Vector_2.h>
#include <CGAL/Cartesian/Direction_2.h>
#include <CGAL/Cartesian/Line_2.h>
#include <CGAL/Cartesian/Ray_2.h>
#include <CGAL/Cartesian/Segment_2.h>
#include <CGAL/Cartesian/Triangle_2.h>
#include <CGAL/Cartesian/Circle_2.h>
#include <CGAL/Cartesian/Iso_rectangle_2.h>
#include <CGAL/Cartesian/Aff_transformation_2.h>
#include <CGAL/Cartesian/Data_accessor_2.h>

#include <CGAL/Cartesian/global_operators_2.h>
#include <CGAL/Cartesian/predicates_on_points_2.h>
#include <CGAL/Cartesian/predicates_on_directions_2.h>
#include <CGAL/Cartesian/predicates_on_lines_2.h>
#include <CGAL/Cartesian/distance_predicates_2.h>

#include <CGAL/Cartesian/constructions_on_points_2.h>
#include <CGAL/Cartesian/constructions_on_lines_2.h>
#include <CGAL/Cartesian/constructions_on_circles_2.h>
#include <CGAL/Cartesian/distance_computations_2.h>

#include <CGAL/Cartesian/Aff_transformation_2.C>

// #include <CGAL/Kernel/Construction_objects_2.h>
// #include <CGAL/Kernel/Predicate_objects_2.h>
#include <CGAL/predicate_objects_on_points_2.h>
#include <CGAL/Kernel/function_objects.h>

CGAL_BEGIN_NAMESPACE

// This class is a restricted 2D geometric kernel
// It is useful only if you do not need the 3D kernel
// If you need both, you should be using Cartesian<FT>

template< class FT_ >
struct Cartesian_2 :
  public Cartesian_base_2< Cartesian_2<FT_>, FT_ >
{
    // Number types and representation tag
    typedef FT_                                 RT;
    typedef FT_                                 FT;
    typedef Cartesian_tag                       Rep_tag;
    typedef Cartesian_tag                       Kernel_tag;

    typedef Cartesian_2<FT_>                    Self;
    typedef Cartesian_base_2<Self,FT_>          Kernel_base;

    typedef typename Kernel_base::Object_2        Object_2;
    
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
    // The classes are inherited and because of partial specialization,
    // Cartesian_2<FT>::Point_2 is exactly CGAL::Point_2< Cartesian_2<FT> >
    // We still need to inherit explicitly, see Cartesian.h for explanation

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

#else
    // Now CGAL::Point_2<R> is only a wrapper around CGAL::PointC2<R>
    // It is necessary to redefine here the classes to ensure that
    // Cartesian_2<FT>::Point_2 is exactly CGAL::Point_2< Cartesian_2<FT> >

    // Cartesian_2<FT>::Base is needed so that CGAL::Point_2< Cartesian_2<FT> >
    // can inherit from Cartesian_2<FT>::Point_2_base

    typedef typename Kernel_base::Point_2                 Point_2_base;
    typedef typename Kernel_base::Vector_2                Vector_2_base;
    typedef typename Kernel_base::Direction_2             Direction_2_base;
    typedef typename Kernel_base::Segment_2               Segment_2_base;
    typedef typename Kernel_base::Line_2                  Line_2_base;
    typedef typename Kernel_base::Ray_2                   Ray_2_base;
    typedef typename Kernel_base::Triangle_2              Triangle_2_base;
    typedef typename Kernel_base::Circle_2                Circle_2_base;
    typedef typename Kernel_base::Iso_rectangle_2         Iso_rectangle_2_base;
    typedef typename Kernel_base::Aff_transformation_2    
                                                     Aff_transformation_2_base;

    // Note: necessary to qualify Point_2 by CGAL:: to disambiguate between
    // Point_2 in the current namespace (nested within CGAL)
    // CGAL::Point_2< Cartesian_2<FT> > (which is in the CGAL namespace)

    typedef CGAL::Point_2<Self>                           Point_2;
    typedef CGAL::Vector_2<Self>                          Vector_2;
    typedef CGAL::Direction_2<Self>                       Direction_2;
    typedef CGAL::Segment_2<Self>                         Segment_2;
    typedef CGAL::Line_2<Self>                            Line_2;
    typedef CGAL::Ray_2<Self>                             Ray_2;
    typedef CGAL::Triangle_2<Self>                        Triangle_2;
    typedef CGAL::Circle_2<Self>                          Circle_2;
    typedef CGAL::Iso_rectangle_2<Self>                   Iso_rectangle_2;
    typedef CGAL::Aff_transformation_2<Self>              Aff_transformation_2;

    typedef Data_accessorC2<Self>                         Data_accessor_2;
    typedef ConicCPA2<Point_2,Data_accessor_2>            Conic_2;

#endif // CGAL_CFG_NO_ADVANCED_KERNEL

    // TODO: cleanup
    static  FT make_FT(const RT & num, const RT& denom) { return num/denom;}
    static  FT make_FT(const RT & num)                  { return num;}
    static  RT FT_numerator(const FT &r)                { return r;}
    static  RT FT_denominator(const FT &)               { return RT(1);}

typedef CGALi::Construct<Point_2>              Construct_point_2;
typedef CGALi::Construct<Vector_2>             Construct_vector_2;
typedef CGALi::Construct<Direction_2>          Construct_direction_2;
typedef CGALi::Construct<Segment_2>            Construct_segment_2;
typedef CGALi::Construct<Line_2>               Construct_line_2;
typedef CGALi::Construct<Ray_2>                Construct_ray_2;
typedef CGALi::Construct<Circle_2>             Construct_circle_2;
typedef CGALi::Construct<Triangle_2>           Construct_triangle_2;
typedef CGALi::Construct<Iso_rectangle_2>      Construct_iso_rectangle_2;
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

Construct_iso_rectangle_2
construct_iso_rectangle_2_object() const
{ return Construct_iso_rectangle_2();}

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

typedef CGALi::Call_direction_to_get<Direction_2> 
                                                Construct_direction_of_line_2;
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

typedef CGALi::Call_perpendicular_to_get<Vector_2> 
                                           Construct_perpendicular_vector_2;
Construct_perpendicular_vector_2
construct_perpendicular_vector_2_object() const 
{ return Construct_perpendicular_vector_2(); }

typedef CGALi::Call_perpendicular_to_get<Direction_2>  
                                       Construct_perpendicular_direction_2;
Construct_perpendicular_direction_2
construct_perpendicular_direction_2_object() const 
{ return Construct_perpendicular_direction_2(); }

typedef CGALi::Call_perpendicular_to_get<Line_2>   
                                        Construct_perpendicular_line_2;
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

typedef CGALi::Call_opposite_to_get<Segment_2>     
                                                Construct_opposite_segment_2;
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

typedef CGALi::Call_opposite_to_get<Triangle_2>    
                                               Construct_opposite_triangle_2;
Construct_opposite_triangle_2
construct_opposite_triangle_2_object() const 
{ return Construct_opposite_triangle_2(); }

typedef CGALi::Call_opposite_to_get<Circle_2>      Construct_opposite_circle_2;
Construct_opposite_circle_2
construct_opposite_circle_2_object() const 
{ return Construct_opposite_circle_2(); }

typedef CGALi::Assign                                  Assign_2;
Assign_2
assign_2_object() const 
{ return Assign_2(); }


typedef CGALi::Call_transform                      Transform_2;
Transform_2
transform_2_object() const 
{ return Transform_2(); }

typedef CGALi::Intersect                           Intersect_2;
Intersect_2
intersect_2_object() const 
{ return Intersect_2(); }

typedef CGALi::Call_y_at_x_to_get<FT>              Compute_y_at_x_2;
Compute_y_at_x_2
compute_y_at_x_2_object() const 
{ return Compute_y_at_x_2(); }

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

typedef CGALi::Equal_xy                            Equal_xy_2;
Equal_xy_2
equal_xy_2_object() const 
{ return Equal_xy_2(); }

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

typedef CGAL ::p_Less_dist_to_point<Point_2>       Less_distance_to_point_2;
Less_distance_to_point_2
less_distance_to_point_2_object(const Point_2& p) const 
{ return Less_distance_to_point_2(p); }

typedef CGAL ::p_Less_dist_to_line_2p<Point_2>  Less_signed_distance_to_line_2;
Less_signed_distance_to_line_2
less_signed_distance_to_line_2_object(const Point_2& p, 
				      const Point_2& q) const 
{ return Less_signed_distance_to_line_2(p,q); }

typedef CGAL ::p_Less_rotate_ccw<Point_2>          Less_rotate_ccw_2;
Less_rotate_ccw_2
less_rotate_ccw_2(const Point_2& p) const 
{ return Less_rotate_ccw_2(p); }

typedef CGALi::Counterclockwise_in_between   Counterclockwise_in_between_2;
Counterclockwise_in_between_2
counterclockwise_in_between_2_object() const 
{ return Counterclockwise_in_between_2(); }

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

typedef CGALi::Are_strictly_ordered_along_line     
                                           Are_strictly_ordered_along_line_2;
Are_strictly_ordered_along_line_2
are_strictly_ordered_along_line_2_object() const 
{ return Are_strictly_ordered_along_line_2(); }


typedef CGALi::Collinear_are_ordered_along_line    
                                            Collinear_are_ordered_along_line_2;
Collinear_are_ordered_along_line_2
collinear_are_ordered_along_line_2_object() const 
{ return Collinear_are_ordered_along_line_2(); }

typedef CGALi::Collinear_are_strictly_ordered_along_line 
                                   Collinear_are_strictly_ordered_along_line_2;
Collinear_are_strictly_ordered_along_line_2
collinear_are_strictly_ordered_along_line_2_object() const 
{ return Collinear_are_strictly_ordered_along_line_2(); }

};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_2_H
