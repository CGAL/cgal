// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Mael Rouxel-Labbé
//                 Konstantinos Katrioplas (konst.katrioplas@gmail.com)

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_MESH_SMOOTHING_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_MESH_SMOOTHING_IMPL_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/internal/Smoothing/smoothing_helpers.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/iterator.h>
#include <CGAL/number_type_config.h>
#include <CGAL/Origin.h>
#include <CGAL/property_map.h>
#include <CGAL/utils.h>

#include <boost/graph/graph_traits.hpp>

#include "ceres/ceres.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <typename V, typename GT>
double get_radian_angle(const V& v1, const V& v2, const GT& gt)
{
  return gt.compute_approximate_angle_3_object()(v1, v2) * CGAL_PI / 180.;
}

// super naive for now. Not sure it even makes sense to do something like that for surfaces
template<typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
class Delaunay_edge_flipper
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor      edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::reference       Point_ref;
  typedef typename GeomTraits::Vector_3                                    Vector;

public:
  Delaunay_edge_flipper(TriangleMesh& mesh,
                        const VertexPointMap vpmap,
                        const GeomTraits& traits)
    : mesh_(mesh), vpmap_(vpmap), traits_(traits)
  { }

  bool should_be_flipped(const edge_descriptor e) const
  {
    if(is_border(e, mesh_))
      return false;

    const halfedge_descriptor h = halfedge(e, mesh_);
    const halfedge_descriptor opp_h = opposite(h, mesh_);

    vertex_descriptor v0 = source(h, mesh_);
    vertex_descriptor v1 = target(h, mesh_);
    vertex_descriptor v2 = target(next(h, mesh_), mesh_);
    vertex_descriptor v3 = target(next(opp_h, mesh_), mesh_);
    const Point_ref p0 = get(vpmap_, v0);
    const Point_ref p1 = get(vpmap_, v1);
    const Point_ref p2 = get(vpmap_, v2);
    const Point_ref p3 = get(vpmap_, v3);

    double alpha = get_radian_angle(Vector(p0 - p2), Vector(p1 - p2), traits_);
    double beta = get_radian_angle(Vector(p1 - p3), Vector(p0 - p3), traits_);

    // not local Delaunay if the sum of the angles is greater than pi
    if(alpha + beta <= CGAL_PI)
      return false;

    // Don't want to flip if the other diagonal already exists
    // @todo remeshing can be used to still flip those
    std::pair<edge_descriptor, bool> other_hd_already_exists = edge(v2, v3, mesh_);
    if(other_hd_already_exists.second)
      return false;

    return true;
  }

  template <typename Marked_edges_map, typename EdgeRange>
  void add_to_stack_if_unmarked(const edge_descriptor e,
                                const Marked_edges_map marks,
                                EdgeRange& edge_range)
  {
    if(!get(marks, e))
    {
      put(marks, e, true);
      edge_range.push_back(e);
    }
  }

public:
  template <typename FaceRange>
  void operator()(const FaceRange& face_range)
  {
    std::cout << " we be flippin'" << std::endl;

    // edges to consider
    std::vector<edge_descriptor> edge_range;
    edge_range.reserve(3 * face_range.size());
    for(face_descriptor f : face_range)
    {
      for(halfedge_descriptor h : halfedges_around_face(halfedge(f, mesh_), mesh_))
        edge_range.push_back(edge(h, mesh_));
    }

    // keep unique elements
    std::sort(edge_range.begin(), edge_range.end());
    edge_range.erase(std::unique(edge_range.begin(), edge_range.end()), edge_range.end());

    std::cout << "unique range of size: " << edge_range.size() << std::endl;

    // Mark edges that are in the stack
    typedef CGAL::dynamic_edge_property_t<bool>                         Edge_property_tag;
    typedef typename boost::property_map<TriangleMesh,
                                         Edge_property_tag>::type       Marked_edges_map;

    Marked_edges_map marks = get(Edge_property_tag(), mesh_);

    // dynamic pmaps do not have default values...
    for(edge_descriptor e : edges(mesh_))
      put(marks, e, false);
    for(edge_descriptor e : edge_range)
      put(marks, e, true);

    int flipped_n = 0;
    while(!edge_range.empty())
    {
      edge_descriptor e = edge_range.back();

      edge_range.pop_back();
      put(marks, e, false);

      if(should_be_flipped(e))
      {
        ++flipped_n;

        halfedge_descriptor h = halfedge(e, mesh_);
        Euler::flip_edge(h, mesh_);

        add_to_stack_if_unmarked(edge(next(h, mesh_), mesh_), marks, edge_range);
        add_to_stack_if_unmarked(edge(prev(h, mesh_), mesh_), marks, edge_range);
        add_to_stack_if_unmarked(edge(next(opposite(h, mesh_), mesh_), mesh_), marks, edge_range);
        add_to_stack_if_unmarked(edge(prev(opposite(h, mesh_), mesh_), mesh_), marks, edge_range);
      }
    }

    std::cout << "Flipped " << flipped_n << " times" << std::endl;
  }

private:
  TriangleMesh& mesh_;
  const VertexPointMap vpmap_;
  const GeomTraits& traits_;
};

template<typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
class Angle_smoother
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::reference       Point_ref;
  typedef typename GeomTraits::Vector_3                                    Vector;

  typedef std::pair<halfedge_descriptor, halfedge_descriptor>              He_pair;

public:
  Angle_smoother(const TriangleMesh& mesh,
                 const VertexPointMap vpmap,
                 const GeomTraits& traits)
    : mesh_(mesh), vpmap_(vpmap), traits_(traits)
  { }

private:
  Vector rotate_edge(const halfedge_descriptor main_he,
                     const He_pair& incident_pair) const
  {
    // get common vertex around which the edge is rotated
    Point_ref pt = get(vpmap_, target(main_he, mesh_));

    Point_ref left_pt = get(vpmap_, source(incident_pair.first, mesh_));
    Point_ref right_pt = get(vpmap_, target(incident_pair.second, mesh_));
    CGAL_assertion(target(incident_pair.first, mesh_) == source(incident_pair.second, mesh_));

    Vector edge1(pt, left_pt);
    Vector edge2(pt, right_pt);

    // find bisector
    internal::normalize(edge1, traits_);
    internal::normalize(edge2, traits_);

    Vector bisector = traits_.construct_sum_of_vectors_3_object()(edge1, edge2);
    internal::normalize(bisector, traits_);

    return bisector;
  }

public:
  // If it's ever allowed to move vertices on the border, the min angle computations will be missing
  // some values (angles incident to the border)
  Vector operator()(const vertex_descriptor v) const
  {
    Vector move = CGAL::NULL_VECTOR;
    double weights_sum = 0.;

    for(halfedge_descriptor main_he : halfedges_around_source(v, mesh_))
    {
      He_pair incident_pair = std::make_pair(prev(opposite(main_he, mesh_), mesh_),
                                             next(main_he, mesh_));

      // avoid zero angles
      Point_ref ps = get(vpmap_, source(main_he, mesh_));
      Point_ref pt = get(vpmap_, target(main_he, mesh_));
      Point_ref left_pt = get(vpmap_, source(incident_pair.first, mesh_));
      Point_ref right_pt = get(vpmap_, target(incident_pair.second, mesh_));
      CGAL_assertion(target(incident_pair.first, mesh_) == source(incident_pair.second, mesh_));

      Vector left_v(pt, left_pt);
      Vector right_v(pt, right_pt);

      // rotate
      double angle = get_radian_angle(right_v, left_v, traits_);
      CGAL_warning(angle != 0.); // no degenerate faces is a precondition
      if(angle == 0.)
        continue;

      Vector bisector = rotate_edge(main_he, incident_pair);
      double scaling_factor = CGAL::approximate_sqrt(
                                traits_.compute_squared_distance_3_object()(get(vpmap_, source(main_he, mesh_)),
                                                                            get(vpmap_, target(main_he, mesh_))));
      bisector = traits_.construct_scaled_vector_3_object()(bisector, scaling_factor);
      Vector ps_psi(ps, traits_.construct_translated_point_3_object()(pt, bisector));

      // small angles carry more weight
      double weight = 1. / (angle*angle);
      weights_sum += weight;

      move += weight * ps_psi;
    }

    if(weights_sum != 0.)
     move /= weights_sum;

    return move;
  }

private:
  const TriangleMesh& mesh_;
  const VertexPointMap vpmap_;
  const GeomTraits& traits_;
};

template<typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
class Area_smoother
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type      Point;
  typedef typename boost::property_traits<VertexPointMap>::reference       Point_ref;
  typedef typename GeomTraits::Vector_3                                    Vector;

public:
  Area_smoother(const TriangleMesh& mesh,
                const VertexPointMap vpmap,
                const GeomTraits& traits)
    : mesh_(mesh), vpmap_(vpmap), traits_(traits)
  { }

private:
  double element_area(const vertex_descriptor v1,
                      const vertex_descriptor v2,
                      const vertex_descriptor v3) const
  {
    return CGAL::to_double(CGAL::approximate_sqrt(traits_.compute_squared_area_3_object()(get(vpmap_, v1),
                                                                                          get(vpmap_, v2),
                                                                                          get(vpmap_, v3))));
  }

  double element_area(const Point& P,
                      const vertex_descriptor v2,
                      const vertex_descriptor v3) const
  {
    return CGAL::to_double(CGAL::approximate_sqrt(traits_.compute_squared_area_3_object()(P,
                                                                                          get(vpmap_, v2),
                                                                                          get(vpmap_, v3))));
  }

  double compute_average_area_around(const vertex_descriptor v) const
  {
    double sum_areas = 0.;
    unsigned int number_of_edges = 0;

    for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
    {
      // opposite vertices
      vertex_descriptor vi = source(next(h, mesh_), mesh_);
      vertex_descriptor vj = target(next(h, mesh_), mesh_);

      double S = element_area(v, vi, vj);
      sum_areas += S;
      ++number_of_edges;
    }

    return sum_areas / number_of_edges;
  }

  struct Face_energy
  {
    Face_energy(const Point& pi, const Point& pj, const double s_av)
      :
        qx(pi.x()), qy(pi.y()), qz(pi.z()),
        rx(pj.x()), ry(pj.y()), rz(pj.z()),
        s_av(s_av)
    { }

    // just for convencience
    template <typename T>
    double area(const T x, const T y, const T z) const {
      return CGAL::approximate_sqrt(CGAL::squared_area(Point(x, y, z),
                                                       Point(qx, qy, qz),
                                                       Point(rx, ry, rz)));
    }

    template <typename T>
    double evaluate(const T x, const T y, const T z) const { return area(x, y, z) - s_av; }

    template <typename T>
    bool operator()(const T* const x, const T* const y, const T* const z,
                    T* residual) const
    {
#define CGAL_CERES_USE_NUMERIC_DIFFERENCIATION
#ifdef CGAL_CERES_USE_NUMERIC_DIFFERENCIATION
      residual[0] = evaluate(x[0], y[0], z[0]);
#else
      // Computations must be explicit so automatic differenciation can be used
      T dqx = qx - x[0];
      T dqy = qy - y[0];
      T dqz = qz - z[0];
      T drx = rx - x[0];
      T dry = ry - y[0];
      T drz = rz - z[0];

      T vx = dqy*drz - dqz*dry;
      T vy = dqz*drx - dqx*drz;
      T vz = dqx*dry - dqy*drx;

      T squared_area = 0.25 * (vx*vx + vy*vy + vz*vz);
      T area = sqrt(squared_area);

      residual[0] = area - s_av;
#endif
      return true;
    }

  private:
    const double qx, qy, qz;
    const double rx, ry, rz;
    const double s_av;
  };

public:
  Vector operator()(const vertex_descriptor v) const
  {
    const Point_ref vp = get(vpmap_, v);

    const double S_av = compute_average_area_around(v);

    const double initial_x = vp.x();
    const double initial_y = vp.y();
    const double initial_z = vp.z();
    double x = initial_x, y = initial_y, z = initial_z;

    ceres::Problem problem;

    for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
    {
      CGAL_assertion(!is_border(h, mesh_));

      vertex_descriptor vi = source(next(h, mesh_), mesh_);
      vertex_descriptor vj = target(next(h, mesh_), mesh_);
      const Point_ref vip = get(vpmap_, vi);
      const Point_ref vjp = get(vpmap_, vj);

#ifdef CGAL_CERES_USE_NUMERIC_DIFFERENCIATION
      ceres::CostFunction* cost_function =
        new ceres::NumericDiffCostFunction<Face_energy, ceres::CENTRAL, 1, 1, 1, 1>(new Face_energy(vip, vjp, S_av));
#else
      ceres::CostFunction* cost_function =
        new ceres::AutoDiffCostFunction<Face_energy, 1, 1, 1, 1>(new Face_energy(vip, vjp, S_av));
#endif
      problem.AddResidualBlock(cost_function, NULL, &x, &y, &z);
    }

    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";
    std::cout << "x : " << initial_x << " -> " << x << "\n";
    std::cout << "y : " << initial_y << " -> " << y << "\n";
    std::cout << "z : " << initial_z << " -> " << z << "\n";

    return Vector(x - initial_x, y - initial_y, z - initial_z);
  }

private:
  const TriangleMesh& mesh_;
  const VertexPointMap vpmap_;
  const GeomTraits& traits_;
};

template<typename Optimizer, typename TriangleMesh,
         typename VertexPointMap, typename VertexConstraintMap,
         typename GeomTraits>
class Mesh_smoother
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor      edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type      Point;
  typedef typename boost::property_traits<VertexPointMap>::reference       Point_ref;
  typedef typename GeomTraits::FT                                          FT;
  typedef typename GeomTraits::Vector_3                                    Vector;
  typedef typename GeomTraits::Segment_3                                   Segment;
  typedef typename GeomTraits::Triangle_3                                  Triangle;

  typedef std::vector<Triangle>                                            Triangle_list;
  typedef std::pair<halfedge_descriptor, halfedge_descriptor>              He_pair;

  typedef std::vector<Triangle>                                            Triangle_container;
  typedef CGAL::AABB_triangle_primitive<GeomTraits, typename Triangle_container::iterator> AABB_Primitive;
  typedef CGAL::AABB_traits<GeomTraits, AABB_Primitive>                    AABB_Traits;
  typedef CGAL::AABB_tree<AABB_Traits>                                     Tree;

public:
  Mesh_smoother(TriangleMesh& pmesh,
                VertexPointMap& vpmap,
                VertexConstraintMap& vcmap,
                const GeomTraits& traits)
    :
      mesh_(pmesh), vpmap_(vpmap), vcmap_(vcmap), traits_(traits)
  {}

  ~Mesh_smoother() { delete tree_ptr_; }

public:
  template<typename FaceRange>
  void set_vertex_range(const FaceRange& face_range)
  {
    vrange_.reserve(3 * face_range.size());
    for(face_descriptor f : face_range)
    {
     for(vertex_descriptor v : vertices_around_face(halfedge(f, mesh_), mesh_))
      vrange_.push_back(v);
    }

    // get rid of duplicate vertices
    std::sort(vrange_.begin(), vrange_.end());
    vrange_.erase(std::unique(vrange_.begin(), vrange_.end()), vrange_.end());
  }

  template<typename FaceRange>
  void init_smoothing(const FaceRange& face_range)
  {
    CGAL_precondition(CGAL::is_triangle_mesh(mesh_));
    CGAL_precondition_code(std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor> degen_faces;)
    CGAL_precondition_code(CGAL::Polygon_mesh_processing::degenerate_faces(
                            mesh_, std::inserter(degen_faces, degen_faces.begin()),
                            CGAL::parameters::vertex_point_map(vpmap_).geom_traits(traits_));)
    CGAL_precondition(degen_faces.empty());

    set_vertex_range(face_range);

    input_triangles_.clear();
    input_triangles_.reserve(face_range.size());

    for(face_descriptor f : face_range)
    {
      halfedge_descriptor h = halfedge(f, mesh_);
      input_triangles_.push_back(traits_.construct_triangle_3_object()(get(vpmap_, source(h, mesh_)),
                                                                       get(vpmap_, target(h, mesh_)),
                                                                       get(vpmap_, target(next(h, mesh_), mesh_))));
    }

    tree_ptr_ = new Tree(input_triangles_.begin(), input_triangles_.end());
    tree_ptr_->accelerate_distance_queries();
  }

  // generic optimizer, the move is computed by 'Optimizer'
  std::size_t optimize(const bool use_sanity_checks = true,
                       const bool apply_moves_in_single_batch = false,
                       const bool enforce_no_min_angle_regression = false)
  {
    typedef CGAL::dynamic_vertex_property_t<Point>                      Vertex_property_tag;
    typedef typename boost::property_map<TriangleMesh,
                                         Vertex_property_tag>::type     Position_map;

    Position_map new_positions = get(Vertex_property_tag(), mesh_);

    Optimizer compute_move(mesh_, vpmap_, traits_);

#ifdef CGAL_PMP_SMOOTHING_DEBUG
    double total_displacement = 0;
#endif

    std::size_t moved_points = 0;
    for(vertex_descriptor v : vrange_)
    {
      if(is_border(v, mesh_) || is_constrained(v))
        continue;

      // compute normal to v
      Vector vn = compute_vertex_normal(v, mesh_, CGAL::parameters::vertex_point_map(vpmap_)
                                                                   .geom_traits(traits_));

      // calculate movement
      const Point_ref pos = get(vpmap_, v);
      Vector move = compute_move(v);

      // @test the effect of this and whether it should be a parameter or not
      if(true/*project_move_on_tangent_plane*/)
      {
        // Gram Schmidt so that the new location is on the tangent plane of v (i.e. do mv -= (mv*n)*n)
        const FT sp = traits_.compute_scalar_product_3_object()(vn, move);
        move = traits_.construct_sum_of_vectors_3_object()(move,
                 traits_.construct_scaled_vector_3_object()(vn, - sp));
      }

      Point new_pos = pos + move;

      // @tmp
//      std::cout << "sanity check? " << use_sanity_checks << " test: " << does_move_create_bad_faces(v, new_pos) << std::endl;
//      std::cout << "min check? " << enforce_no_min_angle_regression << " test: " << does_improve_min_angle_in_star(v, new_pos) << std::endl;

      if((!use_sanity_checks || !does_move_create_bad_faces(v, new_pos)) &&
         (!enforce_no_min_angle_regression || does_improve_min_angle_in_star(v, new_pos)))
      {
#ifdef CGAL_PMP_SMOOTHING_DEBUG
        std::cout << "moving " << get(vpmap_, v) << " to " << new_pos << std::endl;
        total_displacement += CGAL::approximate_sqrt(traits_.compute_squared_length_3_object()(move));
#endif

        if(apply_moves_in_single_batch)
          put(new_positions, v, new_pos);
        else
          put(vpmap_, v, new_pos);

        ++moved_points;
      }
      else // some sanity check failed
      {
#ifdef CGAL_PMP_SMOOTHING_DEBUG
        std::cout << "move is rejected!" << std::endl;
#endif
        if(apply_moves_in_single_batch)
          put(new_positions, v, pos);
      }
    }

    // update locations
    if(apply_moves_in_single_batch)
    {
      for(vertex_descriptor v : vrange_)
      {
        if(is_border(v, mesh_) || is_constrained(v))
          continue;

        put(vpmap_, v, get(new_positions, v));
      }
    }

#ifdef CGAL_PMP_SMOOTHING_DEBUG
    std::cout << "moved: " << moved_points << " points based on angle." << std::endl;
    std::cout << "total displacement: " << total_displacement << std::endl;
    std::cout << "not improved min angle: " << vrange_.size() - moved_points << " times." << std::endl;
#endif

    return moved_points;
  }

  void project_to_surface()
  {
#ifdef CGAL_PMP_SMOOTHING_DEBUG
    std::cout << "Projecting back to the surface" << std::endl;
#endif

    for(vertex_descriptor v : vrange_)
    {
      if(is_border(v, mesh_) || is_constrained(v))
        continue;

      Point_ref p_query = get(vpmap_, v);
      Point projected = tree_ptr_->closest_point(p_query);
      std::cout << p_query << " is projected to: " << projected << std::endl;
      put(vpmap_, v, projected);
    }
  }

private:
  bool is_constrained(const vertex_descriptor v)
  {
    return get(vcmap_, v);
  }

  // check for degenerate or inversed faces
  bool does_move_create_bad_faces(const vertex_descriptor v,
                                  const Point& new_pos) const
  {
    // check for null faces and face inversions
    for(halfedge_descriptor main_he : halfedges_around_source(v, mesh_))
    {
      const halfedge_descriptor prev_he = prev(main_he, mesh_);
      const Point_ref lpt = get(vpmap_, target(main_he, mesh_));
      const Point_ref rpt = get(vpmap_, source(prev_he, mesh_));

      if(traits_.collinear_3_object()(lpt, rpt, new_pos))
        return true;

      const Point_ref old_pos = get(vpmap_, v);
      Vector ov_1 = traits_.construct_vector_3_object()(old_pos, lpt);
      Vector ov_2 = traits_.construct_vector_3_object()(old_pos, rpt);
      Vector old_n = traits_.construct_cross_product_vector_3_object()(ov_1, ov_2);
      Vector nv_1 = traits_.construct_vector_3_object()(new_pos, lpt);
      Vector nv_2 = traits_.construct_vector_3_object()(new_pos, rpt);
      Vector new_n = traits_.construct_cross_product_vector_3_object()(nv_1, nv_2);

      if(!is_positive(traits_.compute_scalar_product_3_object()(old_n, new_n)))
      {
#ifdef CGAL_PMP_SMOOTHING_DEBUG
      std::cout << "Moving vertex would result in the inversion of a face normal!" << std::endl;
#endif
        return true;
      }
    }

    return false;
  }

  bool does_improve_min_angle_in_star(const vertex_descriptor v,
                                      const Point& new_pos) const
  {
    // check if the minimum angle of the star has not deteriorated
    double old_min_angle = CGAL_PI;
    for(halfedge_descriptor main_he : halfedges_around_source(v, mesh_))
    {
      const Point_ref old_pos = get(vpmap_, v);

      const halfedge_descriptor prev_he = prev(main_he, mesh_);
      const Point_ref lpt = get(vpmap_, target(main_he, mesh_));
      const Point_ref rpt = get(vpmap_, source(prev_he, mesh_));

      old_min_angle = (std::min)(old_min_angle,
                                 (std::min)(get_radian_angle(Vector(old_pos, lpt),
                                                             Vector(old_pos, rpt), traits_),
                                            (std::min)(get_radian_angle(Vector(lpt, rpt),
                                                                        Vector(lpt, old_pos), traits_),
                                                       get_radian_angle(Vector(rpt, old_pos),
                                                                        Vector(rpt, lpt), traits_))));
    }

    for(halfedge_descriptor main_he : halfedges_around_source(v, mesh_))
    {
      const halfedge_descriptor prev_he = prev(main_he, mesh_);
      const Point_ref lpt = get(vpmap_, target(main_he, mesh_));
      const Point_ref rpt = get(vpmap_, source(prev_he, mesh_));

      if(get_radian_angle(Vector(new_pos, lpt), Vector(new_pos, rpt), traits_) < old_min_angle ||
         get_radian_angle(Vector(lpt, rpt), Vector(lpt, new_pos), traits_) < old_min_angle ||
         get_radian_angle(Vector(rpt, new_pos), Vector(rpt, lpt), traits_) < old_min_angle)
      {
#ifdef CGAL_PMP_SMOOTHING_DEBUG
        const Point_ref old_pos = get(vpmap_, v);

        std::cout << "deterioration of min angle in the star!" << std::endl;
        std::cout << "old/new positions: " << old_pos << " " << new_pos << std::endl;;
        std::cout << "old min angle: " << old_min_angle << std::endl;
        std::cout << "new angles: " << std::endl;
        std::cout << get_radian_angle(Vector(new_pos, lpt), Vector(new_pos, rpt), traits_) << " ";
        std::cout << get_radian_angle(Vector(lpt, rpt), Vector(lpt, new_pos), traits_) << " ";
        std::cout << get_radian_angle(Vector(rpt, new_pos), Vector(rpt, lpt), traits_) << std::endl;
#endif
        return false;
      }
    }

    return true;
  }

private:
  TriangleMesh& mesh_;
  VertexPointMap& vpmap_;
  VertexConstraintMap vcmap_;
  GeomTraits traits_;

  Tree* tree_ptr_;
  std::vector<vertex_descriptor> vrange_;
  Triangle_container input_triangles_;
};

} // namespace internal
} // namespace Polygon_mesh_processing
} // namespace CGAL


#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_MESH_SMOOTHING_IMPL_H
