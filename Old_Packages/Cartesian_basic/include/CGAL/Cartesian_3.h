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
// file          : include/CGAL/Cartesian_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_3_H
#define CGAL_CARTESIAN_3_H

#include <CGAL/basic.h>
#include <CGAL/cartesian_classes.h>

#ifdef CGAL_CFG_NO_ADVANCED_KERNEL
  // Because we cannot use Michael's scheme, we need the wrapper classes
  // We include them (they are common to Cartesian and Homogeneous)
  #include <CGAL/user_classes.h>
#endif // CGAL_CFG_NO_ADVANCED_KERNEL

#define CGAL_REP_CLASS_DEFINED
#define CGAL_CARTESIAN_CLASS_DEFINED
#define CGAL_NO_2D_IN_3D_KERNEL

CGAL_BEGIN_NAMESPACE

template< class R, class FT_ >
struct Cartesian_base_3
{
    // Number types and representation tag
    typedef FT_                                   RT;
    typedef FT_                                   FT;
    typedef Cartesian_tag                         Rep_tag;
    typedef Cartesian_tag                         Kernel_tag;
    typedef CGAL::Object                          Object_3;
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

    typedef CGAL::Point_3<R,Rep_tag>              Point_3;
    typedef CGAL::Vector_3<R,Rep_tag>             Vector_3;
    typedef CGAL::Direction_3<R,Rep_tag>          Direction_3;
    typedef CGAL::Line_3<R,Rep_tag>               Line_3;
    typedef CGAL::Plane_3<R,Rep_tag>              Plane_3;
    typedef CGAL::Ray_3<R,Rep_tag>                Ray_3;
    typedef CGAL::Segment_3<R,Rep_tag>            Segment_3;
    typedef CGAL::Triangle_3<R,Rep_tag>           Triangle_3;
    typedef CGAL::Tetrahedron_3<R,Rep_tag>        Tetrahedron_3;
    typedef CGAL::Sphere_3<R,Rep_tag>             Sphere_3;
    typedef CGAL::Iso_cuboid_3<R,Rep_tag>         Iso_cuboid_3;
    typedef CGAL::Aff_transformation_3<R,Rep_tag> Aff_transformation_3;
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

    typedef PointC3<R>                            Point_3;
    typedef VectorC3<R>                           Vector_3;
    typedef DirectionC3<R>                        Direction_3;
    typedef LineC3<R>                             Line_3;
    typedef PlaneC3<R>                            Plane_3;
    typedef RayC3<R>                              Ray_3;
    typedef SegmentC3<R>                          Segment_3;
    typedef TriangleC3<R>                         Triangle_3;
    typedef TetrahedronC3<R>                      Tetrahedron_3;
    typedef SphereC3<R>                           Sphere_3;
    typedef Iso_cuboidC3<R>                       Iso_cuboid_3;
    typedef Aff_transformationC3<R>               Aff_transformation_3;
#endif // CGAL_CFG_NO_ADVANCED_KERNEL
};

CGAL_END_NAMESPACE
 
#include <CGAL/Cartesian/Point_2.h>
#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/Cartesian/Vector_3.h>
#include <CGAL/Cartesian/Direction_3.h>

#include <CGAL/Cartesian/Point_2.C>
#include <CGAL/Cartesian/Point_3.C>
#include <CGAL/Cartesian/Vector_3.C>
#include <CGAL/Cartesian/Direction_3.C>

#include <CGAL/Cartesian/Line_3.h>
#include <CGAL/Cartesian/Plane_3.h>
#include <CGAL/Cartesian/Ray_3.h>
#include <CGAL/Cartesian/Segment_3.h>
#include <CGAL/Cartesian/Triangle_3.h>
#include <CGAL/Cartesian/Tetrahedron_3.h>
#include <CGAL/Cartesian/Iso_cuboid_3.h>
#include <CGAL/Cartesian/Aff_transformation_3.h>

#include <CGAL/Cartesian/global_operators_3.h>
#include <CGAL/Cartesian/predicates_on_points_3.h>
#include <CGAL/Cartesian/predicates_on_planes_3.h>
#include <CGAL/Cartesian/distance_predicates_3.h>

#include <CGAL/Cartesian/constructions_on_points_3.h>
#include <CGAL/Cartesian/constructions_on_vectors_3.h>
#include <CGAL/Cartesian/constructions_on_planes_3.h>
#include <CGAL/Cartesian/distance_computations_3.h>

#include <CGAL/Cartesian/Line_3.C>
#include <CGAL/Cartesian/Plane_3.C>
#include <CGAL/Cartesian/Ray_3.C>
#include <CGAL/Cartesian/Segment_3.C>
#include <CGAL/Cartesian/Triangle_3.C>
#include <CGAL/Cartesian/Tetrahedron_3.C>
#include <CGAL/Cartesian/Iso_cuboid_3.C>
#include <CGAL/Cartesian/Aff_transformation_3.C>
 
// #include <CGAL/Kernel/Construction_objects_3.h>
// #include <CGAL/Kernel/Predicate_objects_3.h>
#include <CGAL/predicate_objects_on_points_2.h>
#include <CGAL/Kernel/function_objects.h>

CGAL_BEGIN_NAMESPACE

// This class is a restricted 3D geometric kernel
// It is useful only if you do not need the 2D kernel
// If you need both, you should be using Cartesian<FT>

template< class FT_ >
struct Cartesian_3 :
  public Cartesian_base_3< Cartesian_3<FT_>, FT_ >
{
    // Number types and representation tag
    typedef FT_                                   RT;
    typedef FT_                                   FT;
    typedef Cartesian_tag                         Rep_tag;
    typedef Cartesian_tag                         Kernel_tag;

    typedef Cartesian_3<FT_>                      Self;
    typedef Cartesian_base_3<Self,FT_>            Kernel_base;

    typedef typename Kernel_base::Object_2          Object_2;
    typedef typename Kernel_base::Object_3          Object_3;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
    // The other classes are inherited and because of partial specialization,
    // Cartesian_3<FT>::Point_3 is exactly CGAL::Point_3< Cartesian_3<FT> >
    // We still need to inherit explicitly, see Cartesian.h for explanation

    typedef typename Kernel_base::Point_2         Point_2;
    typedef typename Kernel_base::Vector_2        Vector_2;
    typedef typename Kernel_base::Direction_2     Direction_2;
    typedef typename Kernel_base::Segment_2       Segment_2;
    typedef typename Kernel_base::Line_2          Line_2;
    typedef typename Kernel_base::Ray_2           Ray_2;
    typedef typename Kernel_base::Triangle_2      Triangle_2;
    typedef typename Kernel_base::Circle_2        Circle_2;
    typedef typename Kernel_base::Iso_rectangle_2 Iso_rectangle_2;
    typedef typename Kernel_base::Aff_transformation_2 Aff_transformation_2;

    typedef typename Kernel_base::Point_3         Point_3;
    typedef typename Kernel_base::Vector_3        Vector_3;
    typedef typename Kernel_base::Direction_3     Direction_3;
    typedef typename Kernel_base::Line_3          Line_3;
    typedef typename Kernel_base::Plane_3         Plane_3;
    typedef typename Kernel_base::Ray_3           Ray_3;
    typedef typename Kernel_base::Segment_3       Segment_3;
    typedef typename Kernel_base::Triangle_3      Triangle_3;
    typedef typename Kernel_base::Tetrahedron_3   Tetrahedron_3;
    typedef typename Kernel_base::Sphere_3        Sphere_3;
    typedef typename Kernel_base::Iso_cuboid_3    Iso_cuboid_3;
    typedef typename Kernel_base::Aff_transformation_3 Aff_transformation_3;

 #else
    // Now CGAL::Point_3<R> is only a wrapper around CGAL::PointC3<R>
    // It is necessary to redefine here the classes to ensure that
    // Cartesian_3<FT>::Point_3 is exactly CGAL::Point_3< Cartesian_3<FT> >

    // Cartesian_3<FT>::Base is needed so that CGAL::Point_3< Cartesian_3<FT> >
    // can inherit from Cartesian_3<FT>::typename Kernel_base::Point_3

    typedef typename Kernel_base::Point_2         Point_2_base;
    typedef typename Kernel_base::Vector_2        Vector_2_base;
    typedef typename Kernel_base::Direction_2     Direction_2_base;
    typedef typename Kernel_base::Segment_2       Segment_2_base;
    typedef typename Kernel_base::Line_2          Line_2_base;
    typedef typename Kernel_base::Ray_2           Ray_2_base;
    typedef typename Kernel_base::Triangle_2      Triangle_2_base;
    typedef typename Kernel_base::Circle_2        Circle_2_base;
    typedef typename Kernel_base::Iso_rectangle_2 Iso_rectangle_2_base;
    typedef typename Kernel_base::Aff_transformation_2 
                                                  Aff_transformation_2_base;

    typedef typename Kernel_base::Point_3         Point_3_base;
    typedef typename Kernel_base::Vector_3        Vector_3_base;
    typedef typename Kernel_base::Direction_3     Direction_3_base;
    typedef typename Kernel_base::Line_3          Line_3_base;
    typedef typename Kernel_base::Plane_3         Plane_3_base;
    typedef typename Kernel_base::Ray_3           Ray_3_base;
    typedef typename Kernel_base::Segment_3       Segment_3_base;
    typedef typename Kernel_base::Triangle_3      Triangle_3_base;
    typedef typename Kernel_base::Tetrahedron_3   Tetrahedron_3_base;
    typedef typename Kernel_base::Sphere_3        Sphere_3_base;
    typedef typename Kernel_base::Iso_cuboid_3    Iso_cuboid_3_base;
    typedef typename Kernel_base::Aff_transformation_3    
                                                  Aff_transformation_3_base;

    // Note: necessary to qualify Point_3 by CGAL:: to disambiguate between
    // Point_3 in the current namespace (nested within CGAL) and
    // CGAL::Point_3< Cartesian_3<FT> > (which is in the CGAL namespace)

    typedef CGAL::Point_2<Self>                   Point_2;
    typedef CGAL::Vector_2<Self>                  Vector_2;
    typedef CGAL::Direction_2<Self>               Direction_2;
    typedef CGAL::Segment_2<Self>                 Segment_2;
    typedef CGAL::Line_2<Self>                    Line_2;
    typedef CGAL::Ray_2<Self>                     Ray_2;
    typedef CGAL::Triangle_2<Self>                Triangle_2;
    typedef CGAL::Circle_2<Self>                  Circle_2;
    typedef CGAL::Iso_rectangle_2<Self>           Iso_rectangle_2;
    typedef CGAL::Aff_transformation_2<Self>      Aff_transformation_2;

    typedef Data_accessorC2<Self>                 Data_accessor_2;
    typedef ConicCPA2<Point_2,Data_accessor_2>    Conic_2;

    typedef CGAL::Point_3<Self>                   Point_3;
    typedef CGAL::Vector_3<Self>                  Vector_3;
    typedef CGAL::Direction_3<Self>               Direction_3;
    typedef CGAL::Line_3<Self>                    Line_3;
    typedef CGAL::Plane_3<Self>                   Plane_3;
    typedef CGAL::Ray_3<Self>                     Ray_3;
    typedef CGAL::Segment_3<Self>                 Segment_3;
    typedef CGAL::Triangle_3<Self>                Triangle_3;
    typedef CGAL::Tetrahedron_3<Self>             Tetrahedron_3;
    typedef CGAL::Sphere_3<Self>                  Sphere_3;
    typedef CGAL::Iso_cuboid_3<Self>              Iso_cuboid_3;
    typedef CGAL::Aff_transformation_3<Self>      Aff_transformation_3;

#endif // CGAL_CFG_NO_ADVANCED_KERNEL

    // TODO: cleanup
    static   FT make_FT(const RT & num, const RT& denom) { return num/denom;}
    static   FT make_FT(const RT & num)                  { return num;}
    static   RT FT_numerator(const FT &r)                { return r;}
    static   RT FT_denominator(const FT &)               { return RT(1);}

typedef CGALi::Construct<Point_3>              Construct_point_3;
typedef CGALi::Construct<Vector_3>             Construct_vector_3;
typedef CGALi::Construct<Direction_3>          Construct_direction_3;
typedef CGALi::Construct<Segment_3>            Construct_segment_3;
typedef CGALi::Construct<Plane_3>              Construct_plane_3;
typedef CGALi::Construct<Line_3>               Construct_line_3;
typedef CGALi::Construct<Ray_3>                Construct_ray_3;
typedef CGALi::Construct<Triangle_3>           Construct_triangle_3;
typedef CGALi::Construct<Tetrahedron_3>        Construct_tetrahedron_3;
typedef CGALi::Construct<Iso_cuboid_3>         Construct_iso_cuboid_3;
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

Construct_iso_cuboid_3
construct_iso_cuboid_3_object() const
{return Construct_iso_cuboid_3();}

Construct_aff_transformation_3
construct_aff_transformation_3_object() const 
{ return Construct_aff_transformation_3(); }

typedef CGALi::Call_point_to_get<Point_3>          Construct_point_on_3;
Construct_point_on_3
construct_point_on_3_object() const 
{ return Construct_point_on_3(); }

typedef CGALi::Call_second_point_to_get<Point_3>   Construct_second_point_on_3;
Construct_second_point_on_3
construct_second_point_on_3_object() const 
{ return Construct_second_point_on_3(); }

typedef CGALi::Call_perpendicular_plane_to_get<Plane_3> 
                                              Construct_perpendicular_plane_3;
Construct_perpendicular_plane_3
construct_perpendicular_plane_3() const 
{ return Construct_perpendicular_plane_3(); }

typedef CGALi::p_Midpoint<Point_3>                 Construct_midpoint_3;
Construct_midpoint_3
construct_midpoint_3_object() const 
{ return Construct_midpoint_3(); }

typedef CGALi::p_Circumcenter<Point_3>             Construct_circumcenter_3;
Construct_circumcenter_3
construct_circumcenter_3_object() const 
{ return Construct_circumcenter_3(); }

typedef CGALi::Call_opposite_to_get<Segment_3>     
                                              Construct_opposite_segment_3;
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

typedef CGALi::Assign                                  Assign_3;
Assign_3
assign_3_object() const 
{ return Assign_3(); }


typedef CGALi::Intersect                           Intersect_3;
Intersect_3
intersect_3_object() const 
{ return Intersect_3(); }

typedef CGALi::Call_squared_length_to_get<FT>      Compute_squared_length_3;
Compute_squared_length_3
compute_squared_length_3_object() const 
{ return Compute_squared_length_3(); }

typedef CGALi::Equal                               Equal_3;
Equal_3
equal_3_object() const 
{ return Equal_3(); }

typedef CGALi::Equal_x                             Equal_x_3;
Equal_x_3
equal_x_3_object() const 
{ return Equal_x_3(); }

typedef CGALi::Equal_y                             Equal_y_3;
Equal_y_3
equal_y_3_object() const 
{ return Equal_y_3(); }

typedef CGALi::Equal_z                             Equal_z_3;
Equal_z_3
equal_z_3_object() const 
{ return Equal_z_3(); }

typedef CGALi::Equal_xy                            Equal_xy_3;
Equal_xy_3
equal_xy_3_object() const 
{ return Equal_xy_3(); }

typedef CGALi::Equal_xyz                           Equal_xyz_3;
Equal_xyz_3
equal_xyz_3_object() const 
{ return Equal_xyz_3(); }

typedef CGALi::Less_x                              Less_x_3;
Less_x_3
less_x_3_object() const 
{ return Less_x_3(); }

typedef CGALi::Less_y                              Less_y_3;
Less_y_3
less_y_3_object() const 
{ return Less_y_3(); }

typedef CGALi::Less_z                              Less_z_3;
Less_z_3
less_z_3_object() const 
{ return Less_z_3(); }

typedef CGAL::p_Less_xy<Point_3>                   Less_xy_3;
Less_xy_3
less_xy_3_object() const 
{ return Less_xy_3(); }

typedef CGALi::Less_xyz                            Less_xyz_3;
Less_xyz_3
less_xyz_3_object() const 
{ return Less_xyz_3(); }

typedef CGALi::Compare_x                           Compare_x_3;
Compare_x_3
compare_x_3_object() const 
{ return Compare_x_3(); }

typedef CGALi::Compare_y                           Compare_y_3;
Compare_y_3
compare_y_3_object() const 
{ return Compare_y_3(); }

typedef CGALi::Compare_z                           Compare_z_3;
Compare_z_3
compare_z_3_object() const 
{ return Compare_z_3(); }

typedef CGALi::Compare_xy                          Compare_xy_3;
Compare_xy_3
compare_xy_3_object() const 
{ return Compare_xyz_3(); }

typedef CGALi::Compare_xyz                         Compare_xyz_3;
Compare_xyz_3
compare_xyz_3_object() const 
{ return Compare_xyz_3(); }

typedef CGAL ::p_Less_dist_to_point<Point_3>       Less_distance_to_point_3;
Less_distance_to_point_3
less_distance_to_point_3_object(const Point_3& p) const 
{ return Less_distance_to_point_3(p); }

typedef CGALi::Collinear                           Collinear_3;
Collinear_3
collinear_3_object() const 
{ return Collinear_3(); }

typedef CGALi::Coplanar                            Coplanar_3 ;
Coplanar_3
coplanar_3_object() const 
{ return Coplanar_3(); }

typedef CGAL ::p_Orientation<Point_3>              Orientation_3;
Orientation_3
orientation_3_object() const 
{ return Orientation_3(); }

typedef CGALi::Call_is_degenerate                  Is_degenerate_3;
Is_degenerate_3
is_degenerate_3_object() const 
{ return Is_degenerate_3(); }

typedef CGALi::Call_has_on                         Has_on_3;
Has_on_3
has_on_3_object() const 
{ return Has_on_3(); }

typedef CGALi::Call_has_on_bounded_side            Has_on_bounded_side_3;
Has_on_bounded_side_3
has_on_bounded_side_3_object() const 
{ return Has_on_bounded_side_3(); }

typedef CGALi::Call_has_on_unbounded_side          Has_on_unbounded_side_3;
Has_on_unbounded_side_3
has_on_unbounded_side_3_object() const 
{ return Has_on_unbounded_side_3(); }

typedef CGALi::Call_has_on_boundary                Has_on_boundary_3;
Has_on_boundary_3
has_on_boundary_3_object() const 
{ return Has_on_boundary_3(); }

typedef CGALi::Call_has_on_positive_side           Has_on_positive_side_3;
Has_on_positive_side_3
has_on_positive_side_3_object() const 
{ return Has_on_positive_side_3(); }

typedef CGALi::Call_has_on_negative_side           Has_on_negative_side_3;
Has_on_negative_side_3
has_on_negative_side_3_object() const 
{ return Has_on_negative_side_3(); }

typedef CGALi::Call_oriented_side                  Oriented_side_3;
Oriented_side_3
oriented_side_3_object() const 
{ return Oriented_side_3(); }

typedef CGALi::Are_ordered_along_line              Are_ordered_along_line_3;
Are_ordered_along_line_3
are_ordered_along_line_3_object() const 
{ return Are_ordered_along_line_3(); }

typedef CGALi::Are_strictly_ordered_along_line     
                                           Are_strictly_ordered_along_line_3;
Are_strictly_ordered_along_line_3
are_strictly_ordered_along_line_3_object() const 
{ return Are_strictly_ordered_along_line_3(); }

typedef CGALi::Collinear_are_ordered_along_line    
                                          Collinear_are_ordered_along_line_3;
Collinear_are_ordered_along_line_3
collinear_are_ordered_along_line_3_object() const 
{ return Collinear_are_ordered_along_line_3(); }

typedef CGALi::Collinear_are_strictly_ordered_along_line 
                                 Collinear_are_strictly_ordered_along_line_3;
Collinear_are_strictly_ordered_along_line_3
collinear_are_strictly_ordered_along_line_3_object() const 
{ return Collinear_are_strictly_ordered_along_line_3(); }

typedef CGALi::Side_of_oriented_sphere             Side_of_oriented_sphere_3;
Side_of_oriented_sphere_3
side_of_oriented_sphere_3_object() const 
{ return Side_of_oriented_sphere_3(); }

typedef CGALi::Side_of_bounded_sphere              Side_of_bounded_sphere_3;
Side_of_bounded_sphere_3
side_of_bounded_sphere_3_object() const 
{ return Side_of_bounded_sphere_3(); }

};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_3_H
