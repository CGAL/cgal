// Copyright (c) 2000,2001
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel

#ifndef CGAL_CARTESIAN_D_H
#define CGAL_CARTESIAN_D_H

#include <CGAL/basic.h>
#include <CGAL/representation_tags.h>
#include <CGAL/Kernel_d/function_objects.h>
#include <CGAL/Linear_algebraCd.h>
#include <vector>

#include <CGAL/Kernel_d/Kernel_classesCd.h>
#include <CGAL/Kernel_d/PointCd.h>
#include <CGAL/Kernel_d/VectorCd.h>
#include <CGAL/Kernel_d/DirectionCd.h>
#include <CGAL/Kernel_d/HyperplaneCd.h>
#include <CGAL/Kernel_d/Aff_transformationCd.h>
#include <CGAL/Kernel_d/PointCd_impl.h>
#include <CGAL/Kernel_d/VectorCd_impl.h>
#include <CGAL/Kernel_d/DirectionCd_impl.h>
#include <CGAL/Kernel_d/HyperplaneCd_impl.h>
#include <CGAL/Kernel_d/function_objectsCd.h>
#include <CGAL/Kernel_d/intersection_objectsCd.h>
#include <CGAL/Kernel_d/Interface_classes.h>
#include <CGAL/Kernel_d/simple_objects.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template <class pFT, class pLA = Linear_algebraCd<pFT> >
class Cartesian_d
{
public:
  typedef Cartesian_d<pFT,pLA> Self;
  typedef pFT                    RT;
  typedef pFT                    FT;
  typedef pLA                    LA;

  typedef Cartesian_tag        Rep_tag;
  typedef Cartesian_tag        Kernel_tag;

  enum { Has_filtered_predicates = false };
  typedef Boolean_tag<Has_filtered_predicates> Has_filtered_predicates_tag;

  typedef PointCd2<RT,LA>             Point_d_base;
  // renamed because of clash with Cartesian...
  typedef VectorCd<RT,LA>             Vector_d_base;
  typedef DirectionCd<RT,LA>          Direction_d_base;
  typedef HyperplaneCd<RT,LA>         Hyperplane_d_base;
  typedef Aff_transformationCd<RT,LA> Aff_transformation_d_base;

  typedef CGAL::Point_d<Self>              Point_d;
  typedef CGAL::Vector_d<Self>             Vector_d;
  typedef CGAL::Direction_d<Self>          Direction_d;
  typedef CGAL::Hyperplane_d<Self>         Hyperplane_d;
  typedef CGAL::Aff_transformation_d<Self> Aff_transformation_d;

  typedef typename Point_d_base::Cartesian_const_iterator Cartesian_const_iterator_d;

    // Boolean   had originally been Bool. It was renamed to avoid a conflict
    // between a macro defined in Xlib.h poorly chosen to have the same name,
    // that is 'Bool'.
    typedef typename Same_uncertainty_nt<bool, FT>::type
                                                        Boolean;
    typedef typename Same_uncertainty_nt<CGAL::Sign, FT>::type
                                                        Sign;
    typedef typename Same_uncertainty_nt<CGAL::Comparison_result, FT>::type
                                                        Comparison_result;
    typedef typename Same_uncertainty_nt<CGAL::Orientation, FT>::type
                                                        Orientation;
    typedef typename Same_uncertainty_nt<CGAL::Oriented_side, FT>::type
                                                        Oriented_side;
    typedef typename Same_uncertainty_nt<CGAL::Bounded_side, FT>::type
                                                        Bounded_side;
    typedef typename Same_uncertainty_nt<CGAL::Angle, FT>::type
                                                        Angle;

  typedef Dynamic_dimension_tag            Dimension;

  template <typename T>
  struct Ambient_dimension {
    typedef typename T::Ambient_dimension type;
  };

  template <typename T>
  struct Feature_dimension {
    typedef typename T::Feature_dimension type;
  };

  template <typename K>
  class Construct_cartesian_const_iterator
  {
    typedef typename K::Point_d Point_d;
    typedef typename K::Cartesian_const_iterator_d  Cartesian_const_iterator_d;

  public:
    typedef Cartesian_const_iterator_d result_type;

    Cartesian_const_iterator_d
    operator()(const Point_d& p) const
    {
      return p.cartesian_begin();
    }

    Cartesian_const_iterator_d
    operator()( const Point_d& p, int) const
    {
      return p.cartesian_end();
    }
  };

  // TODO: Make it work for the other values
 template <typename K>
  class Construct_vertex
  {
    typedef typename K::Point_d Point_d;
    typedef typename K::Iso_box_d Iso_box_d;
    typedef typename K::Cartesian_const_iterator_d  Cartesian_const_iterator_d;
  public:
    typedef Point_d result_type;

    Point_d operator()(const Iso_box_d&  b, int i)
    {
      if(i == 0){
        return (b.min)();
      }
      return (b.max)();
    }
  };


  typedef Construct_vertex<Self> Construct_vertex_d;

    typedef Construct_cartesian_const_iterator<Self>
                           Construct_cartesian_const_iterator_d;

  Construct_cartesian_const_iterator_d
  construct_cartesian_const_iterator_d_object() const
  {
    return Construct_cartesian_const_iterator_d();
  }

 template <typename K>
  class Construct_min_vertex
  {
    typedef typename K::Point_d Point_d;
    typedef typename K::Iso_box_d Iso_box_d;
  public:
    typedef Point_d result_type;

    Point_d operator()(const Iso_box_d&  b)
    {
      return (b.min)();
    }
  };
  typedef Construct_min_vertex<Self> Construct_min_vertex_d;

  Construct_min_vertex_d
  construct_min_vertex_d_object() const
  {
    return Construct_min_vertex_d();
  }


 template <typename K>
  class Construct_max_vertex
  {
    typedef typename K::Point_d Point_d;
    typedef typename K::Iso_box_d Iso_box_d;
  public:
    typedef Point_d result_type;

    Point_d operator()(const Iso_box_d&  b)
    {
      return (b.max)();
    }
  };
  typedef Construct_max_vertex<Self> Construct_max_vertex_d;

  Construct_max_vertex_d
  construct_max_vertex_d_object() const
  {
    return Construct_max_vertex_d();
  }

  // meta types (fit both kernels):
  typedef CGAL::Sphere_d<Self>   Sphere_d;
  typedef CGAL::Iso_box_d<Self>  Iso_box_d;
  typedef CGAL::Segment_d<Self>  Segment_d;
  typedef CGAL::Ray_d<Self>      Ray_d;
  typedef CGAL::Line_d<Self>     Line_d;

  // construction objects:
  typedef internal::Construct<Point_d> Construct_point_d;
  Construct_point_d construct_point_d_object() const
  { return Construct_point_d(); }

  typedef internal::Construct<Vector_d> Construct_vector_d;
  Construct_vector_d construct_vector_d_object() const
  { return Construct_vector_d(); }

  typedef internal::Construct<Direction_d> Construct_direction_d;
  Construct_direction_d construct_direction_d_object() const
  { return Construct_direction_d(); }

  typedef internal::Construct<Segment_d> Construct_segment_d;
  Construct_segment_d construct_segment_d_object() const
  { return Construct_segment_d(); }

  typedef internal::Construct<Ray_d> Construct_ray_d;
  Construct_ray_d construct_ray_d_object() const
  { return Construct_ray_d(); }

  typedef internal::Construct<Line_d> Construct_line_d;
  Construct_line_d construct_line_d_object() const
  { return Construct_line_d(); }

  typedef internal::Construct<Iso_box_d> Construct_iso_box_d;
  Construct_iso_box_d construct_iso_box_d_object() const
  { return Construct_iso_box_d(); }

  typedef internal::Construct<Hyperplane_d> Construct_hyperplane_d;
  Construct_hyperplane_d construct_hyperplane_d_object() const
  { return Construct_hyperplane_d(); }

  typedef internal::Construct<Sphere_d> Construct_sphere_d;
  Construct_sphere_d construct_sphere_d_object() const
  { return Construct_sphere_d(); }

  typedef internal::Construct<Aff_transformation_d>
    Construct_aff_transformation_d;
  Construct_aff_transformation_d
    construct_aff_transformation_d_object() const
  { return Construct_aff_transformation_d(); }

  // function objects:
  typedef Compute_coordinateCd<Self> Compute_coordinate_d;
  typedef Lift_to_paraboloidCd<Self> Lift_to_paraboloid_d;
  typedef Project_along_d_axisCd<Self> Project_along_d_axis_d;
  typedef MidpointCd<Self> Midpoint_d;
  typedef Squared_distanceCd<Self> Squared_distance_d;
  typedef Position_on_lineCd<Self> Position_on_line_d;
  typedef Barycentric_coordinatesCd<Self> Barycentric_coordinates_d;
  typedef OrientationCd<Self> Orientation_d;
  typedef Coaffine_orientationCd<Self> Coaffine_orientation_d;
  typedef Side_of_oriented_sphereCd<Self> Side_of_oriented_sphere_d;
  typedef Side_of_oriented_subsphereCd<Self> Side_of_oriented_subsphere_d;
  typedef Side_of_bounded_sphereCd<Self> Side_of_bounded_sphere_d;
  typedef Contained_in_simplexCd<Self> Contained_in_simplex_d;
  typedef Contained_in_affine_hullCd<Self> Contained_in_affine_hull_d;
  typedef Affine_rankCd<Self> Affine_rank_d;
  typedef Affinely_independentCd<Self> Affinely_independent_d;
  typedef Compare_lexicographicallyCd<Self> Compare_lexicographically_d;
  typedef Lt_from_compare<Self> Less_lexicographically_d;
  typedef Le_from_compare<Self> Less_or_equal_lexicographically_d;
  typedef Less_coordinateCd<Self> Less_coordinate_d;
  typedef Point_dimensionCd<Self> Point_dimension_d;
  typedef Eq_from_method<Self> Equal_d;
  typedef Center_of_sphereCd<Self> Center_of_sphere_d;
  typedef Contained_in_linear_hullCd<Self> Contained_in_linear_hull_d;
  typedef Linear_rankCd<Self> Linear_rank_d;
  typedef Linearly_independentCd<Self> Linearly_independent_d;
  typedef Linear_baseCd<Self> Linear_base_d;

  Compute_coordinate_d compute_coordinate_d_object() const
  { return Compute_coordinate_d(); }
  Point_dimension_d point_dimension_d_object() const
  { return Point_dimension_d(); }
  Less_coordinate_d less_coordinate_d_object() const
  { return Less_coordinate_d(); }
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
  Coaffine_orientation_d coaffine_orientation_d_object() const
  { return Coaffine_orientation_d(); }
  Side_of_oriented_sphere_d side_of_oriented_sphere_d_object() const
  { return Side_of_oriented_sphere_d(); }
  Side_of_oriented_subsphere_d side_of_oriented_subsphere_d_object() const
  { return Side_of_oriented_subsphere_d(); }
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
  typedef CGAL::Line_line_intersectionCd<Self> Line_line_intersection_d;
  typedef CGAL::Line_hyperplane_intersectionCd<Self>
                                               Line_hyperplane_intersection_d;
  typedef CGAL::Line_d_Line_d_pair<Self> Line_d_Line_d_pair;
  typedef CGAL::Ray_d_Ray_d_pair<Self> Ray_d_Ray_d_pair;
  typedef CGAL::Segment_d_Segment_d_pair<Self> Segment_d_Segment_d_pair;
  typedef CGAL::Line_d_Ray_d_pair<Self> Line_d_Ray_d_pair;
  typedef CGAL::Line_d_Segment_d_pair<Self> Line_d_Segment_d_pair;
  typedef CGAL::Ray_d_Segment_d_pair<Self> Ray_d_Segment_d_pair;
  typedef CGAL::Line_d_Hyperplane_d_pair<Self> Line_d_Hyperplane_d_pair;
  typedef CGAL::Ray_d_Hyperplane_d_pair<Self> Ray_d_Hyperplane_d_pair;
  typedef CGAL::Segment_d_Hyperplane_d_pair<Self> Segment_d_Hyperplane_d_pair;

  typedef internal::Intersect<Self> Intersect_d;
  Intersect_d intersect_d_object() const
  { return Intersect_d(); }

  typedef internal::Do_intersect<Self> Do_intersect_d;
  Do_intersect_d do_intersect_d_object() const
  { return Do_intersect_d(); }

  // FT - RT conversion and access :

  static  FT make_FT(const RT & num, const RT& denom)
  { return num/denom; }

  static  FT make_FT(const RT & num)
  { return FT(num); }

  static  RT FT_numerator(const FT &r)
  { return r; }

  static  RT FT_denominator(const FT & /*r*/)
  { return RT(1); }

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

  typedef internal::Call_has_on_positive_side Has_on_positive_side_d;
  Has_on_positive_side_d has_on_positive_side_d_object() const
  { return Has_on_positive_side_d(); }

  typedef internal::Call_oriented_side Oriented_side_d;
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

}; // Cartesian_d<R>


} //namespace CGAL

#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Kernel_d/Vector_d.h>
#include <CGAL/Kernel_d/Direction_d.h>
#include <CGAL/Kernel_d/Hyperplane_d.h>
#include <CGAL/Kernel_d/Aff_transformation_d.h>
#include <CGAL/Kernel_d/Sphere_d.h>
#include <CGAL/Kernel_d/Iso_box_d.h>
#include <CGAL/Kernel_d/Segment_d.h>
#include <CGAL/Kernel_d/Ray_d.h>
#include <CGAL/Kernel_d/Line_d.h>
#include <CGAL/Kernel_d/Line_d_impl.h>
#include <CGAL/intersections_d.h>

#endif // CGAL_CARTESIAN_D_H
