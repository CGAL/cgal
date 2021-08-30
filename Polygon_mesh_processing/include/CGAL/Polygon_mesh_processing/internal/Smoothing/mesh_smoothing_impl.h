// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé
//                 Konstantinos Katrioplas (konst.katrioplas@gmail.com)

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_MESH_SMOOTHING_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_MESH_SMOOTHING_IMPL_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#ifdef CGAL_PMP_USE_CERES_SOLVER
#include <CGAL/Polygon_mesh_processing/internal/Smoothing/ceres_support.h>
#endif

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>

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
#include <CGAL/use.h>
#include <CGAL/utils.h>

#include <boost/graph/graph_traits.hpp>

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
typename GT::FT get_radian_angle(const V& v1, const V& v2, const GT& gt)
{
  typedef typename GT::FT FT;

  return gt.compute_approximate_angle_3_object()(v1, v2) * CGAL_PI / FT(180);
}

// super naive for now. Not sure it even makes sense to do something like that for surfaces
template<typename TriangleMesh,
         typename VertexPointMap,
         typename ECMap,
         typename GeomTraits>
class Delaunay_edge_flipper
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor      edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::reference       Point_ref;
  typedef typename GeomTraits::FT                                          FT;
  typedef typename GeomTraits::Vector_3                                    Vector;

public:
  Delaunay_edge_flipper(TriangleMesh& mesh,
                        const VertexPointMap vpmap,
                        const ECMap ecmap,
                        const GeomTraits& traits)
    : mesh_(mesh), vpmap_(vpmap), ecmap_(ecmap), traits_(traits)
  { }

  bool should_be_flipped(const edge_descriptor e) const
  {
    if(is_border(e, mesh_) || get(ecmap_, e))
      return false;

    const halfedge_descriptor h = halfedge(e, mesh_);
    const halfedge_descriptor opp_h = opposite(h, mesh_);

    const vertex_descriptor v0 = source(h, mesh_);
    const vertex_descriptor v1 = target(h, mesh_);
    const vertex_descriptor v2 = target(next(h, mesh_), mesh_);
    const vertex_descriptor v3 = target(next(opp_h, mesh_), mesh_);

    std::set<vertex_descriptor> unique_vs { v0, v1, v2, v3 };
    if(unique_vs.size() != 4)
      return false;

    // Don't want to flip if the other diagonal already exists
    // @todo remeshing can be used to still flip those
    std::pair<edge_descriptor, bool> other_hd_already_exists = edge(v2, v3, mesh_);
    if(other_hd_already_exists.second)
      return false;

    // not local Delaunay := sum of the opposite angles is greater than pi
    const Point_ref p0 = get(vpmap_, v0);
    const Point_ref p1 = get(vpmap_, v1);
    const Point_ref p2 = get(vpmap_, v2);
    const Point_ref p3 = get(vpmap_, v3);

    FT alpha = get_radian_angle(Vector(p0 - p2), Vector(p1 - p2), traits_);
    FT beta = get_radian_angle(Vector(p1 - p3), Vector(p0 - p3), traits_);

    return (alpha + beta > CGAL_PI);
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
#ifdef CGAL_PMP_SMOOTHING_DEBUG
    std::cout << "Flipping edges" << std::endl;
#endif

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

#ifdef CGAL_PMP_SMOOTHING_DEBUG_PP
        std::cout << "Flipping " << edge(h, mesh_) << std::endl;
#endif
        Euler::flip_edge(h, mesh_);

        add_to_stack_if_unmarked(edge(next(h, mesh_), mesh_), marks, edge_range);
        add_to_stack_if_unmarked(edge(prev(h, mesh_), mesh_), marks, edge_range);
        add_to_stack_if_unmarked(edge(next(opposite(h, mesh_), mesh_), mesh_), marks, edge_range);
        add_to_stack_if_unmarked(edge(prev(opposite(h, mesh_), mesh_), mesh_), marks, edge_range);
      }
    }

#ifdef CGAL_PMP_SMOOTHING_DEBUG
    std::cout << flipped_n << " flips" << std::endl;
#endif
  }

private:
  TriangleMesh& mesh_;
  const VertexPointMap vpmap_;
  const ECMap ecmap_;
  const GeomTraits& traits_;
};

template<typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
class Angle_smoother
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::reference       Point_ref;
  typedef typename GeomTraits::FT                                          FT;
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
    FT weights_sum = FT(0);

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
      Vector bisector = rotate_edge(main_he, incident_pair);
      FT scaling_factor = CGAL::approximate_sqrt(
                            traits_.compute_squared_distance_3_object()(get(vpmap_, source(main_he, mesh_)),
                                                                        get(vpmap_, target(main_he, mesh_))));
      bisector = traits_.construct_scaled_vector_3_object()(bisector, scaling_factor);
      Vector ps_psi(ps, traits_.construct_translated_point_3_object()(pt, bisector));

      FT angle = get_radian_angle(right_v, left_v, traits_);
      if(angle == FT(0))
      {
        // no degenerate faces is a precondition, angle can be 0 but it should be a numerical error
        CGAL_warning(!is_degenerate_triangle_face(face(main_he, mesh_), mesh_));

        return ps_psi; // since a small angle gives more weight, a null angle give priority (?)
      }

      // small angles carry more weight
      FT weight = 1. / CGAL::square(angle);
      weights_sum += weight;

      move += weight * ps_psi;
    }

    if(weights_sum != FT(0))
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
  typedef typename GeomTraits::FT                                          FT;
  typedef typename GeomTraits::Vector_3                                    Vector;

public:
  Area_smoother(const TriangleMesh& mesh,
                const VertexPointMap vpmap,
                const GeomTraits& traits)
    : mesh_(mesh), vpmap_(vpmap), traits_(traits)
  { }

private:
  FT element_area(const vertex_descriptor v1,
                  const vertex_descriptor v2,
                  const vertex_descriptor v3) const
  {
    return CGAL::approximate_sqrt(traits_.compute_squared_area_3_object()(get(vpmap_, v1),
                                                                          get(vpmap_, v2),
                                                                          get(vpmap_, v3)));
  }

  FT element_area(const Point& P,
                  const vertex_descriptor v2,
                  const vertex_descriptor v3) const
  {
    return CGAL::approximate_sqrt(traits_.compute_squared_area_3_object()(P,
                                                                          get(vpmap_, v2),
                                                                          get(vpmap_, v3)));
  }

  FT compute_average_area_around(const vertex_descriptor v) const
  {
    FT sum_areas = 0;
    unsigned int number_of_edges = 0;

    for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
    {
      // opposite vertices
      vertex_descriptor vi = source(next(h, mesh_), mesh_);
      vertex_descriptor vj = target(next(h, mesh_), mesh_);

      FT S = element_area(v, vi, vj);
      sum_areas += S;
      ++number_of_edges;
    }

    return sum_areas / number_of_edges;
  }

  struct Face_energy
  {
    Face_energy(const Point& pi, const Point& pj, const FT s_av)
      :
        qx(pi.x()), qy(pi.y()), qz(pi.z()),
        rx(pj.x()), ry(pj.y()), rz(pj.z()),
        s_av(s_av)
    { }

    // next two functions are just for convencience, the only thing ceres cares about is the operator()
    template <typename T>
    FT area(const T x, const T y, const T z) const
    {
      return CGAL::approximate_sqrt(CGAL::squared_area(Point(x, y, z),
                                                       Point(qx, qy, qz),
                                                       Point(rx, ry, rz)));
    }

    template <typename T>
    FT evaluate(const T x, const T y, const T z) const { return area(x, y, z) - s_av; }

    template <typename T>
    bool operator()(const T* const x, const T* const y, const T* const z,
                    T* residual) const
    {
      // Defining this because I haven't found much difference empirically (auto-diff being maybe
      // a couple % faster), but numeric differenciation should be stronger in the face
      // of difficult cases. Leaving the auto-differenciation formulation in case somebody really
      // cares about the extra speed.
#define CGAL_CERES_USE_NUMERIC_DIFFERENCIATION

#ifdef CGAL_CERES_USE_NUMERIC_DIFFERENCIATION
      residual[0] = evaluate(x[0], y[0], z[0]);
#else
      // Computations must be explicit so that automatic differenciation can be used
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
    const FT qx, qy, qz;
    const FT rx, ry, rz;
    const FT s_av;
  };

public:
  Vector operator()(const vertex_descriptor v) const
  {
#ifdef CGAL_PMP_USE_CERES_SOLVER
    const Point_ref vp = get(vpmap_, v);

    const FT S_av = compute_average_area_around(v);

    const FT initial_x = vp.x();
    const FT initial_y = vp.y();
    const FT initial_z = vp.z();
    FT x = initial_x, y = initial_y, z = initial_z;

    ceres::Problem problem;

    problem.AddParameterBlock(&x, 1);
    problem.AddParameterBlock(&y, 1);
    problem.AddParameterBlock(&z, 1);

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
    options.logging_type = ceres::SILENT;
//    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
//    std::cout << summary.BriefReport() << "\n";
//    std::cout << "x : " << initial_x << " -> " << x << "\n";
//    std::cout << "y : " << initial_y << " -> " << y << "\n";
//    std::cout << "z : " << initial_z << " -> " << z << "\n";

    return Vector(x - initial_x, y - initial_y, z - initial_z);
#else
    CGAL_USE(v);
    return CGAL::NULL_VECTOR;
#endif
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

public:
  Mesh_smoother(TriangleMesh& pmesh,
                VertexPointMap& vpmap,
                VertexConstraintMap& vcmap,
                const GeomTraits& traits)
    :
      mesh_(pmesh), vpmap_(vpmap), vcmap_(vcmap), traits_(traits)
  {}

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
    FT total_displacement = 0;
    std::cout << "apply_moves_in_single_batch: " << apply_moves_in_single_batch << std::endl;
#endif

    std::size_t moved_points = 0;
    for(vertex_descriptor v : vrange_)
    {
      if(is_border(v, mesh_) || is_constrained(v))
        continue;

#ifdef CGAL_PMP_SMOOTHING_DEBUG_PP
      std::cout << "Considering " << v << " pos: " << get(vpmap_, v) << std::endl;
#endif

      // compute normal to v
      Vector vn = compute_vertex_normal(v, mesh_, CGAL::parameters::vertex_point_map(vpmap_)
                                                                   .geom_traits(traits_));

      // calculate movement
      const Point_ref pos = get(vpmap_, v);
      Vector move = compute_move(v);

      // Gram Schmidt so that the new location is on the tangent plane of v (i.e. do mv -= (mv*n)*n)
      const FT sp = traits_.compute_scalar_product_3_object()(vn, move);
      move = traits_.construct_sum_of_vectors_3_object()(
               move, traits_.construct_scaled_vector_3_object()(vn, - sp));

      const Point new_pos = pos + move;
      if(move != CGAL::NULL_VECTOR &&
         !does_move_create_degenerate_faces(v, new_pos) &&
         (!use_sanity_checks || !does_move_create_bad_faces(v, new_pos)) &&
         (!enforce_no_min_angle_regression || does_improve_min_angle_in_star(v, new_pos)))
      {
#ifdef CGAL_PMP_SMOOTHING_DEBUG_PP
        std::cout << "moving " << get(vpmap_, v) << " to " << new_pos << std::endl;
#endif

        if(apply_moves_in_single_batch)
          put(new_positions, v, new_pos);
        else
          put(vpmap_, v, new_pos);

#ifdef CGAL_PMP_SMOOTHING_DEBUG
        total_displacement += CGAL::approximate_sqrt(traits_.compute_squared_length_3_object()(move));
#endif

        ++moved_points;
      }
      else // some sanity check failed
      {
#ifdef CGAL_PMP_SMOOTHING_DEBUG_PP
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

  template <typename AABBTree>
  void project_to_surface(const AABBTree& tree)
  {
#ifdef CGAL_PMP_SMOOTHING_DEBUG
    std::cout << "Projecting back to the surface" << std::endl;
#endif

    for(vertex_descriptor v : vrange_)
    {
      if(is_border(v, mesh_) || is_constrained(v))
        continue;

      Point_ref p_query = get(vpmap_, v);
      const Point projected = tree.closest_point(p_query);
#ifdef CGAL_PMP_SMOOTHING_DEBUG_PP
      std::cout << p_query << " to " << projected << std::endl;
#endif

      put(vpmap_, v, projected);
    }
  }

private:
  bool is_constrained(const vertex_descriptor v)
  {
    return get(vcmap_, v);
  }

  // Null faces are bad because they make normal computation difficult
  bool does_move_create_degenerate_faces(const vertex_descriptor v,
                                         const Point& new_pos) const
  {
    for(halfedge_descriptor main_he : halfedges_around_source(v, mesh_))
    {
      const halfedge_descriptor prev_he = prev(main_he, mesh_);
      const Point_ref lpt = get(vpmap_, target(main_he, mesh_));
      const Point_ref rpt = get(vpmap_, source(prev_he, mesh_));

      if(traits_.collinear_3_object()(lpt, rpt, new_pos))
        return true;
    }

    return false;
  }

  // check for degenerate or inversed faces
  bool does_move_create_bad_faces(const vertex_descriptor v,
                                  const Point& new_pos) const
  {
    // check for face inversions
    for(halfedge_descriptor main_he : halfedges_around_source(v, mesh_))
    {
      const halfedge_descriptor prev_he = prev(main_he, mesh_);
      const Point_ref lpt = get(vpmap_, target(main_he, mesh_));
      const Point_ref rpt = get(vpmap_, source(prev_he, mesh_));

      CGAL_assertion(!traits_.collinear_3_object()(lpt, rpt, new_pos)); // checked above

      const Point_ref old_pos = get(vpmap_, v);
      Vector ov_1 = traits_.construct_vector_3_object()(old_pos, lpt);
      Vector ov_2 = traits_.construct_vector_3_object()(old_pos, rpt);
      Vector old_n = traits_.construct_cross_product_vector_3_object()(ov_1, ov_2);
      Vector nv_1 = traits_.construct_vector_3_object()(new_pos, lpt);
      Vector nv_2 = traits_.construct_vector_3_object()(new_pos, rpt);
      Vector new_n = traits_.construct_cross_product_vector_3_object()(nv_1, nv_2);

      if(!is_positive(traits_.compute_scalar_product_3_object()(old_n, new_n)))
      {
#ifdef CGAL_PMP_SMOOTHING_DEBUG_PP
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
    FT old_min_angle = CGAL_PI;
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
#ifdef CGAL_PMP_SMOOTHING_DEBUG_PP
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

  std::vector<vertex_descriptor> vrange_;
};

} // namespace internal
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_MESH_SMOOTHING_IMPL_H
