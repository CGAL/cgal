// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_MESH_3_POLYHEDRAL_TO_LABELED_FUNCTION_WRAPPER_H
#define CGAL_MESH_3_POLYHEDRAL_TO_LABELED_FUNCTION_WRAPPER_H

#include <CGAL/license/Mesh_3.h>



#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/AABB_polyhedron_segment_primitive.h>

#include <CGAL/point_generators_3.h>

namespace CGAL {

namespace Mesh_3 {


namespace details {

inline
double
max_length(const Bbox_3& b)
{
  return (std::max)(b.xmax()-b.xmin(),
                    (std::max)(b.ymax()-b.ymin(),b.zmax()-b.zmin()) );
}

template <class K>
inline
typename K::Sphere_3
bounding_sphere(const Bbox_3& bbox)
{
  const double radius = max_length(bbox)*2;

  typename K::Point_3 center( (bbox.xmin()+bbox.xmax())/2.,
                              (bbox.ymin()+bbox.ymax())/2.,
                              (bbox.zmin()+bbox.zmax())/2. );

  return typename K::Sphere_3(center, radius * radius + 1.);
}


} // end namespace details

template<class Polyhedron_, class BGT>
class Polyhedral_to_labeled_function_wrapper
{
public:
  // AABB_tree Types
  typedef class AABB_const_polyhedron_triangle_primitive<BGT, Polyhedron_>
                                                    AABB_primitive;
  typedef class AABB_traits<BGT,AABB_primitive>     AABB_traits;
  typedef class AABB_tree<AABB_traits>              AABB_tree;
  typedef typename AABB_traits::Bounding_box        Bounding_box;

  // Geometric traits types
  typedef typename BGT::Point_3     Point_3;
  typedef typename BGT::Sphere_3    Sphere_3;

  // Return type
  typedef int return_type;

  /// Constructor
  Polyhedral_to_labeled_function_wrapper(const Polyhedron_& p,
                                         const int nb_step,
                                         const double step_size,
                                         const double first_level)
    : tree_(new AABB_tree(p.facets_begin(), p.facets_end()))
    , hint_(p.facets_begin()->halfedge()->vertex()->point())
    , nb_step_(nb_step)
    , step_size_(step_size)
    , first_level_(first_level)
  {

  }

  // Default copy constructor and assignment operator are ok

  /// Destructor
  ~Polyhedral_to_labeled_function_wrapper() {}

  /// Operator ()
  return_type operator()(const Point_3& p, bool use_cache=false) const
  {
    const double sq_distance = sq_sign_dist(p,use_cache);

    double lvl = first_level_;
    double sq_lvl = CGAL::sign(lvl)*lvl*lvl;
    for ( int i = 1 ; i <= nb_step_ ; ++i )
    {
      if ( sq_distance < sq_lvl ) return i;
      lvl += step_size_;
      sq_lvl = CGAL::sign(lvl)*lvl*lvl;
    }
    return 0;
  }

  Sphere_3 bounding_sphere() const
  {
    return details::bounding_sphere<BGT>(tree_->bbox());
  }

private:
  double sq_sign_dist(const Point_3& p, bool use_cache) const
  {
    const Bounding_box bbox = tree_->bbox();
    const double diameter = details::max_length(bbox) * 2;

    typedef typename BGT::RT RT;
    typedef typename BGT::FT FT;
    typedef typename BGT::Segment_3 Segment_3;

    Random_points_on_sphere_3<Point_3> random_point(1.);

    typename BGT::Construct_vector_3 vector =
      BGT().construct_vector_3_object();
    typename BGT::Construct_segment_3 segment =
      BGT().construct_segment_3_object();
    typename BGT::Construct_translated_point_3 translate =
      BGT().construct_translated_point_3_object();
    typename BGT::Construct_scaled_vector_3 scale =
      BGT().construct_scaled_vector_3_object();

    Segment_3 seg =
      segment(p,
              translate(p, scale(vector(ORIGIN,*random_point), diameter)) );

    double sign = 1.;
    if (   p.x() >= bbox.xmin() && p.x() <= bbox.xmax()
        && p.y() >= bbox.ymin() && p.y() <= bbox.ymax()
        && p.z() >= bbox.zmin() && p.z() <= bbox.zmax()
        && ((tree_->number_of_intersected_primitives(seg) % 2) == 1) )
    {
        // Point is inside domain
        sign = -1;
    }

    Point_3 proj;
    if ( use_cache )
      proj = tree_->closest_point(p,hint_);
    else
      proj = tree_->closest_point(p);

    // update cache
    hint_ = proj;

    // Return sign squared distance
    typename BGT::Compute_squared_distance_3 distance =
      BGT().compute_squared_distance_3_object();

    return (sign*CGAL::to_double(distance(p,proj)));
  }

private:
  /// Functions to wrap
  AABB_tree* tree_;
  mutable Point_3 hint_;
  const int nb_step_;
  const double step_size_;
  const double first_level_;

};  // end class Polyhedral_to_labeled_function_wrapper



template<class Polyhedron_, class BGT>
class Polyhedral_tolerance_to_labeled_function_wrapper
{
public:
  // AABB_tree Types
  typedef class AABB_const_polyhedron_triangle_primitive<BGT, Polyhedron_>
                                                    AABB_primitive;
  typedef class AABB_traits<BGT,AABB_primitive>     AABB_traits;
  typedef class AABB_tree<AABB_traits>              AABB_tree;
  typedef typename AABB_traits::Bounding_box        Bounding_box;

  // Geometric traits types
  typedef typename BGT::Point_3     Point_3;
  typedef typename BGT::Sphere_3    Sphere_3;

  // Return type
  typedef int return_type;

  /// Constructor
  Polyhedral_tolerance_to_labeled_function_wrapper(const Polyhedron_& p,
                                                   const double tolerance_size)
    : tree_(new AABB_tree(p.facets_begin(), p.facets_end()))
    , sq_tolerance_size_(tolerance_size*tolerance_size)
    , hint_(p.facets_begin()->halfedge()->vertex()->point())
  {
  }

  // Default copy constructor and assignment operator are ok

  /// Destructor
  ~Polyhedral_tolerance_to_labeled_function_wrapper() {}

  /// Operator ()
  return_type operator()(const Point_3& p, bool use_cache=false) const
  {
    if ( sq_dist(p,use_cache) < sq_tolerance_size_ )
      return 1;
    else
      return 0;
  }

  Sphere_3 bounding_sphere() const
  {
    return details::bounding_sphere<BGT>(tree_->bbox());
  }

private:
  double sq_dist(const Point_3& p, bool use_cache) const
  {
    Point_3 proj;
    if ( use_cache )
      proj = tree_->closest_point(p,hint_);
    else
      proj = tree_->closest_point(p);

    // update cache
    hint_ = proj;

    // Return squared distance
    typename BGT::Compute_squared_distance_3 distance =
      BGT().compute_squared_distance_3_object();

    return (CGAL::to_double(distance(p,proj)));
  }

private:
  /// Functions to wrap
  AABB_tree* tree_;
  const double sq_tolerance_size_;
  mutable Point_3 hint_;

};  // end class Polyhedral_tolerance_to_labeled_function_wrapper



template<class Polyhedron_, class BGT>
class Polyhedral_edge_tolerance_to_labeled_function_wrapper
{
public:
  // AABB_tree Types
  typedef class AABB_const_polyhedron_edge_primitive<BGT, Polyhedron_>
                                                    AABB_primitive;
  typedef class AABB_traits<BGT,AABB_primitive>     AABB_traits;
  typedef class AABB_tree<AABB_traits>              AABB_tree;
  typedef typename AABB_traits::Bounding_box        Bounding_box;

  // Geometric traits types
  typedef typename BGT::Point_3     Point_3;
  typedef typename BGT::Sphere_3    Sphere_3;

  // Return type
  typedef int return_type;

  /// Constructor
  Polyhedral_edge_tolerance_to_labeled_function_wrapper(const Polyhedron_& p,
                                                        const double tolerance_size)
    : tree_(new AABB_tree(p.edges_begin(), p.edges_end()))
    , sq_tolerance_size_(tolerance_size*tolerance_size)
    , hint_(p.facets_begin()->halfedge()->vertex()->point())
  {
  }

  // Default copy constructor and assignment operator are ok

  /// Destructor
  ~Polyhedral_edge_tolerance_to_labeled_function_wrapper() {}

  /// Operator ()
  return_type operator()(const Point_3& p, bool use_cache=false) const
  {
    if ( sq_dist(p,use_cache) < sq_tolerance_size_ )
      return 1;
    else
      return 0;
  }

  Sphere_3 bounding_sphere() const
  {
    return details::bounding_sphere<BGT>(tree_->bbox());
  }

private:
  double sq_dist(const Point_3& p, bool use_cache) const
  {
    Point_3 proj;
    if ( use_cache )
      proj = tree_->closest_point(p,hint_);
    else
      proj = tree_->closest_point(p);

    // update cache
    hint_ = proj;

    // Return squared distance
    typename BGT::Compute_squared_distance_3 distance =
      BGT().compute_squared_distance_3_object();

    return (CGAL::to_double(distance(p,proj)));
  }

private:
  /// Functions to wrap
  AABB_tree* tree_;
  const double sq_tolerance_size_;
  mutable Point_3 hint_;

};  // end class Polyhedral_edge_tolerance_to_labeled_function_wrapper

}  // end namespace Mesh_3

}  // end namespace CGAL

#endif // CGAL_MESH_3_POLYHEDRAL_TO_LABELED_FUNCTION_WRAPPER_H
