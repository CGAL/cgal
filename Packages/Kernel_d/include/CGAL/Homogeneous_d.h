// ======================================================================
//
// Copyright (c) 2000,2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Homogeneous_d.h
// package       : Kernel_d
// maintainer    : Michael Seel <Michael.Seel@mpi-sb.mpg.de>
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Seel
// coordinator   : MPI Saarbruecken (Susan.Hert@mpi-sb.mpg.de)
//
// ======================================================================

#ifndef HOMOGENEOUS_D_H
#define HOMOGENEOUS_D_H

#ifndef NOCGALINCL
#include <CGAL/basic.h>
#include <CGAL/Quotient.h>
#include <CGAL/representation_tags.h>
#include <CGAL/Kernel/function_objects.h>
#include <CGAL/Linear_algebraHd.h>
#include <vector>
#endif
#include <CGAL/Kernel_d/Kernel_classesHd.h>
#include <CGAL/Kernel_d/PointHd.h>
#include <CGAL/Kernel_d/VectorHd.h>
#include <CGAL/Kernel_d/DirectionHd.h>
#include <CGAL/Kernel_d/HyperplaneHd.h>
#include <CGAL/Kernel_d/Aff_transformationHd.h>
#include <CGAL/Kernel_d/PointHd.C>
#include <CGAL/Kernel_d/VectorHd.C>
#include <CGAL/Kernel_d/DirectionHd.C>
#include <CGAL/Kernel_d/HyperplaneHd.C>
#include <CGAL/Kernel_d/function_objectsHd.h>
#include <CGAL/Kernel_d/intersection_objectsHd.h>
#include <CGAL/Kernel_d/Interface_classes.h>
#include <CGAL/Kernel_d/simple_objects.h>

CGAL_BEGIN_NAMESPACE

template <class pRT, class pLA = Linear_algebraHd<pRT> >
class Homogeneous_d 
{
public:
  typedef Homogeneous_d<pRT,pLA> Self;
  typedef pRT                    RT;
  typedef Quotient<RT>           FT;
  typedef pLA                    LA;

  typedef Homogeneous_tag        Rep_tag;

  typedef PointHd<RT,LA>              Point_d_base;
  typedef VectorHd<RT,LA>             Vector_d_base;
  typedef DirectionHd<RT,LA>          Direction_d_base;
  typedef HyperplaneHd<RT,LA>         Hyperplane_d_base;
  typedef Aff_transformationHd<RT,LA> Aff_transformation_d_base;

  typedef CGAL::Point_d<Self>              Point_d;
  typedef CGAL::Vector_d<Self>             Vector_d;
  typedef CGAL::Direction_d<Self>          Direction_d;
  typedef CGAL::Hyperplane_d<Self>         Hyperplane_d;
  typedef CGAL::Aff_transformation_d<Self> Aff_transformation_d;

  // meta types (fit both kernels):
  typedef CGAL::Sphere_d<Self>   Sphere_d;
  typedef CGAL::Segment_d<Self>  Segment_d;
  typedef CGAL::Ray_d<Self>      Ray_d;
  typedef CGAL::Line_d<Self>     Line_d;

  // construction objects:
  typedef CGALi::Construct<Point_d> Construct_point_d;
  Construct_point_d construct_point_d_object() const
  { return Construct_point_d(); }

  typedef CGALi::Construct<Vector_d> Construct_vector_d;
  Construct_vector_d construct_vector_d_object() const
  { return Construct_vector_d(); }

  typedef CGALi::Construct<Direction_d> Construct_direction_d;
  Construct_direction_d construct_direction_d_object() const
  { return Construct_direction_d(); }

  typedef CGALi::Construct<Segment_d> Construct_segment_d;
  Construct_segment_d construct_segment_d_object() const
  { return Construct_segment_d(); }

  typedef CGALi::Construct<Ray_d> Construct_ray_d;
  Construct_ray_d construct_ray_d_object() const
  { return Construct_ray_d(); }

  typedef CGALi::Construct<Line_d> Construct_line_d;
  Construct_line_d construct_line_d_object() const
  { return Construct_line_d(); }

  typedef CGALi::Construct<Hyperplane_d> Construct_hyperplane_d;
  Construct_hyperplane_d construct_hyperplane_d_object() const
  { return Construct_hyperplane_d(); }

  typedef CGALi::Construct<Sphere_d> Construct_sphere_d;
  Construct_sphere_d construct_sphere_d_object() const
  { return Construct_sphere_d(); }

  typedef CGALi::Construct<Aff_transformation_d> 
    Construct_aff_transformation_d;
  Construct_aff_transformation_d 
    construct_aff_transformation_d_object() const
  { return Construct_aff_transformation_d(); }

  // function objects:
  typedef Lift_to_paraboloidHd<Self> Lift_to_paraboloid_d;
  typedef Project_along_d_axisHd<Self> Project_along_d_axis_d;
  typedef MidpointHd<Self> Midpoint_d;
  typedef Squared_distanceHd<Self> Squared_distance_d;
  typedef Position_on_lineHd<Self> Position_on_line_d;
  typedef Barycentric_coordinatesHd<Self> Barycentric_coordinates_d;
  typedef OrientationHd<Self> Orientation_d;
  typedef Side_of_oriented_sphereHd<Self> Side_of_oriented_sphere_d;
  typedef Side_of_bounded_sphereHd<Self> Side_of_bounded_sphere_d;
  typedef Contained_in_simplexHd<Self> Contained_in_simplex_d;
  typedef Contained_in_affine_hullHd<Self> Contained_in_affine_hull_d;
  typedef Affine_rankHd<Self> Affine_rank_d;
  typedef Affinely_independentHd<Self> Affinely_independent_d;
  typedef Compare_lexicographicallyHd<Self> Compare_lexicographically_d;
  typedef Lt_from_compare<Self> Less_lexicographically_d;
  typedef Le_from_compare<Self> Less_or_equal_lexicographically_d;
  typedef Eq_from_method<Self> Equal_d;
  typedef Center_of_sphereHd<Self> Center_of_sphere_d;
  typedef Contained_in_linear_hullHd<Self> Contained_in_linear_hull_d;
  typedef Linear_rankHd<Self> Linear_rank_d;
  typedef Linearly_independentHd<Self> Linearly_independent_d;
  typedef Linear_baseHd<Self> Linear_base_d;

  Lift_to_paraboloid_d lift_to_paraboloid_d_object() const
  { return Lift_to_paraboloid_d(); }
  Project_along_d_axis_d project_along_d_axis_d_object() const
  { return Project_along_d_axis_d(); }
  Midpoint_d midpoint_d_object() const
  { return Midpoint_d(); }
  Squared_distance_d squared_distance_d_object() const
  { return Squared_distance_d(); }
  Position_on_line_d position_on_line_d_object() const
  { return Position_on_line_d(); }
  Barycentric_coordinates_d barycentric_coordinates_d_object() const
  { return Barycentric_coordinates_d(); }
  Orientation_d orientation_d_object() const
  { return Orientation_d(); }
  Side_of_oriented_sphere_d side_of_oriented_sphere_d_object() const
  { return Side_of_oriented_sphere_d(); }
  Side_of_bounded_sphere_d side_of_bounded_sphere_d_object() const
  { return Side_of_bounded_sphere_d(); }
  Contained_in_simplex_d contained_in_simplex_d_object() const
  { return Contained_in_simplex_d(); }
  Contained_in_affine_hull_d contained_in_affine_hull_d_object() const
  { return Contained_in_affine_hull_d(); }
  Affine_rank_d affine_rank_d_object() const 
  { return Affine_rank_d(); }
  Affinely_independent_d affinely_independent_d_object() const
  { return Affinely_independent_d(); }
  Equal_d equal_d_object() const
  { return Equal_d(); }
  Compare_lexicographically_d compare_lexicographically_d_object() const
  { return Compare_lexicographically_d(); }
  Less_lexicographically_d less_lexicographically_d_object() const
  { return Less_lexicographically_d(); }
  Less_or_equal_lexicographically_d 
    less_or_equal_lexicographically_d_object() const
  { return Less_or_equal_lexicographically_d(); }
  Center_of_sphere_d center_of_sphere_d_object() const
  { return Center_of_sphere_d(); }
  Contained_in_linear_hull_d contained_in_linear_hull_d_object() const
  { return Contained_in_linear_hull_d(); }
  Linear_rank_d linear_rank_d_object() const 
  { return Linear_rank_d(); }
  Linearly_independent_d linearly_independent_d_object() const
  { return Linearly_independent_d(); }
  Linear_base_d linear_base_d_object() const 
  { return Linear_base_d(); }

  // Intersection objects:
  typedef Line_line_intersectionHd<Self> Line_line_intersection_d;
  typedef Line_hyperplane_intersectionHd<Self> Line_hyperplane_intersection_d;
  typedef Line_d_Line_d_pair<Self> Line_d_Line_d_pair;
  typedef Ray_d_Ray_d_pair<Self> Ray_d_Ray_d_pair;
  typedef Segment_d_Segment_d_pair<Self> Segment_d_Segment_d_pair;
  typedef Line_d_Ray_d_pair<Self> Line_d_Ray_d_pair;
  typedef Line_d_Segment_d_pair<Self> Line_d_Segment_d_pair;
  typedef Ray_d_Segment_d_pair<Self> Ray_d_Segment_d_pair;
  typedef Line_d_Hyperplane_d_pair<Self> Line_d_Hyperplane_d_pair;
  typedef Ray_d_Hyperplane_d_pair<Self> Ray_d_Hyperplane_d_pair;
  typedef Segment_d_Hyperplane_d_pair<Self> Segment_d_Hyperplane_d_pair;

  typedef CGALi::Intersect Intersect_d;
  Intersect_d intersect_d_object() const 
  { return Intersect_d(); }

  // FT - RT conversion and access :

  static  FT  make_FT(const RT & num, const RT& denom)
  { return FT(num, denom); }
  
  static  FT  make_FT(const RT & num)
  { return FT(num); }
  
  static  RT FT_numerator(const FT &r)
  { return r.numerator(); }
  
  static  RT FT_denominator(const FT &r)
  { return r.denominator(); }

  // special stuff for traits class character :

  struct Component_accessor_d {
    template <typename C>
    int dimension(const C& c) const { return c.dimension(); }
    template <typename C>
    RT homogeneous(const C& c, int i) { return c.homogeneous(i); }
    template <typename C>
    FT cartesian(const C& c, int i) { return c.cartesian(i); }
  };
  Component_accessor_d component_accessor_d_object() const
  { return Component_accessor_d(); }

  typedef CGALi::Call_has_on_positive_side Has_on_positive_side_d;
  Has_on_positive_side_d has_on_positive_side_d_object() const
  { return Has_on_positive_side_d(); }

  typedef CGALi::Call_oriented_side Oriented_side_d;
  Oriented_side_d oriented_side_d_object() const
  { return Oriented_side_d(); }

  struct Value_at_d {
    RT operator()(const Hyperplane_d& h, const Point_d& p) const
    { return h.value_at(p); }
  };
  Value_at_d value_at_d_object() const
  { return Value_at_d(); }

  struct Point_to_vector_d {
    Vector_d operator()(const Point_d& p) const
    { return p-CGAL::ORIGIN; }
  };
  Point_to_vector_d point_to_vector_d_object() const
  { return Point_to_vector_d(); }

  struct Vector_to_point_d {
    Point_d operator()(const Vector_d& v) const
    { return CGAL::ORIGIN+v; }
  };
  Vector_to_point_d vector_to_point_d_object() const
  { return Vector_to_point_d(); }
  
  struct Orthogonal_vector_d {
    Vector_d operator()(const Hyperplane_d& h) const
    { return h.orthogonal_vector(); }
  };
  Orthogonal_vector_d orthogonal_vector_d_object() const
  { return Orthogonal_vector_d(); }

  struct Point_of_sphere_d {
    Point_d operator()(const Sphere_d& S, int i)
    { return S.point(i); }
  };
  Point_of_sphere_d point_of_sphere_d_object() const
  { return Point_of_sphere_d(); }

}; // Homogeneous_d<R>
 
 
CGAL_END_NAMESPACE

#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Kernel_d/Vector_d.h>
#include <CGAL/Kernel_d/Direction_d.h>
#include <CGAL/Kernel_d/Hyperplane_d.h>
#include <CGAL/Kernel_d/Aff_transformation_d.h>
#include <CGAL/Kernel_d/Sphere_d.h>
#include <CGAL/Kernel_d/Segment_d.h>
#include <CGAL/Kernel_d/Ray_d.h>
#include <CGAL/Kernel_d/Line_d.h>
#include <CGAL/Kernel_d/Line_d.C>
#include <CGAL/intersections_d.h>

#endif // HOMOGENEOUS_D_H

