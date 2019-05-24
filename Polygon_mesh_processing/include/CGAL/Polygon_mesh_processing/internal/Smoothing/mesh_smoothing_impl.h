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
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/internal/Smoothing/smoothing_helpers.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/iterator.h>
#include <CGAL/property_map.h>
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

template<typename PolygonMesh, typename VertexPointMap, typename VertexConstraintMap, typename GeomTraits>
class Compatible_smoother
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor      edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor      face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type     Point_3;
  typedef typename boost::property_traits<VertexPointMap>::reference      Point_ref;
  typedef typename GeomTraits::FT                                         FT;
  typedef typename GeomTraits::Vector_3                                   Vector;
  typedef typename GeomTraits::Segment_3                                  Segment;
  typedef typename GeomTraits::Triangle_3                                 Triangle;

  typedef std::vector<Triangle>                                           Triangle_list;
  typedef std::pair<halfedge_descriptor, halfedge_descriptor>             He_pair;

  typedef CGAL::AABB_face_graph_triangle_primitive<PolygonMesh,
                                                   VertexPointMap,
                                                   CGAL::Tag_true/*single mesh*/,
                                                   CGAL::Tag_true/*cache data*/>
                                                                          AABB_Primitive;
  typedef CGAL::AABB_traits<GeomTraits, AABB_Primitive>                   AABB_Traits;
  typedef CGAL::AABB_tree<AABB_Traits>                                    Tree;

public:
  Compatible_smoother(PolygonMesh& pmesh,
                      VertexPointMap& vpmap,
                      VertexConstraintMap& vcmap,
                      const GeomTraits& traits)
    :
      mesh_(pmesh), vpmap_(vpmap), vcmap_(vcmap), traits_(traits)
  {}

  ~Compatible_smoother() { delete tree_ptr_; }

  template<typename FaceRange>
  void init_smoothing(const FaceRange& face_range)
  {
    set_vertex_range(face_range);

    tree_ptr_ = new Tree(faces(mesh_).begin(), faces(mesh_).end(), mesh_, vpmap_);
    tree_ptr_->accelerate_distance_queries();
  }

  void angle_relaxation()
  {
    typedef CGAL::dynamic_vertex_property_t<Point_3>                    Vertex_property_tag;
    typedef typename boost::property_map<PolygonMesh,
                                         Vertex_property_tag>::type     Position_map;

    Position_map new_positions = get(Vertex_property_tag(), mesh_);

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

      // Gram Schmidt so that the new location is on the tangent plane of v (i.e. mv -= (mv*n)*n)
      const FT sp = traits_.compute_scalar_product_3_object()(vn, move);
      move = traits_.construct_sum_of_vectors_3_object()(move,
               traits_.construct_scaled_vector_3_object()(vn, - sp));

      Point_3 new_pos = pos + move;
      if(does_improve(v, new_pos))
      {
        ++moved_points;
        put(new_positions, v, new_pos);
      }
      else
      {
        std::cout << "move rejected!!" << std::endl;
        put(new_positions, v, pos);
      }
    }

    std::cout << moved_points << " moves" << std::endl;

    // update locations
    for(vertex_descriptor v : vrange_)
    {
      if(is_border(v, mesh_) || is_constrained(v))
        continue;

      put(vpmap_, v, get(new_positions, v));
    }

#ifdef CGAL_PMP_SMOOTHING_DEBUG
    std::cout << "moved: " << moved_points << " points based on angle." << std::endl;
    std::cout << "not improved min angle: " << vrange_.size() - moved_points << " times." << std::endl;
#endif
  }

  void area_relaxation(const double precision)
  {
    std::size_t moved_points = 0;
    for(vertex_descriptor v : vrange_)
    {
       if(!is_border(v, mesh_) && !is_constrained(v))
       {
         if(gradient_descent(v, precision))
           ++moved_points;
       }
    }

#ifdef CGAL_PMP_SMOOTHING_DEBUG
    std::cout << "moved : " << moved_points << " points based on area." << std::endl;
    std::cout << "non convex energy found: " << vrange_.size() - moved_points << " times." << std::endl;
#endif
  }

  void project_to_surface()
  {
    for(vertex_descriptor v : vrange_)
    {
      if(!is_border(v, mesh_) && !is_constrained(v))
      {
        Point_ref p_query = get(vpmap_, v);
        Point_3 projected = tree_ptr_->closest_point(p_query);
        put(vpmap_, v, projected);
      }
    }
  }

private:
  // angle bisecting functions
  // -------------------------
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

  // If it's ever allowed to move vertices on the border, the min angle computations will be missing
  // some values (angles incident to the border)
  Vector compute_move(const vertex_descriptor v) const
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
      double angle = get_angle(right_v, left_v);
      CGAL_assertion(angle > 0.); // no degenerate faces is a precondition

      Vector bisector = rotate_edge(main_he, incident_pair);
      bisector = traits_.construct_scaled_vector_3_object()(bisector,
                                                            CGAL::approximate_sqrt(sqlength(main_he, mesh_)));
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

  // angle measurement & evaluation
  // ------------------------------

  double get_angle(const Vector& e1, const Vector& e2) const
  {
    return traits_.compute_approximate_angle_3_object()(e1, e2) * CGAL_PI / 180.;
  }

  bool does_improve(const vertex_descriptor v,
                    const Point_3& new_pos) const
  {
    std::cout << "DOES IMPROVE AT V " << v << std::endl;

    // check for null faces and face inversions
    for(halfedge_descriptor main_he : halfedges_around_source(v, mesh_))
    {
      const Point_ref old_pos = get(vpmap_, v);

      const halfedge_descriptor prev_he = prev(main_he, mesh_);
      const Point_ref lpt = get(vpmap_, target(main_he, mesh_));
      const Point_ref rpt = get(vpmap_, source(prev_he, mesh_));

      if(traits_.collinear_3_object()(lpt, rpt, new_pos))
        return false;

      Vector ov_1 = traits_.construct_vector_3_object()(old_pos, lpt);
      Vector ov_2 = traits_.construct_vector_3_object()(old_pos, rpt);
      Vector old_n = traits_.construct_cross_product_vector_3_object()(ov_1, ov_2);
      Vector nv_1 = traits_.construct_vector_3_object()(new_pos, lpt);
      Vector nv_2 = traits_.construct_vector_3_object()(new_pos, rpt);
      Vector new_n = traits_.construct_cross_product_vector_3_object()(nv_1, nv_2);

      if(!is_positive(traits_.compute_scalar_product_3_object()(old_n, new_n)))
        return false;
    }

    // check if the minimum angle of the star has not deteriorated
    double old_min_angle = CGAL_PI;
    for(halfedge_descriptor main_he : halfedges_around_source(v, mesh_))
    {
      const Point_ref old_pos = get(vpmap_, v);

      const halfedge_descriptor prev_he = prev(main_he, mesh_);
      const Point_ref lpt = get(vpmap_, target(main_he, mesh_));
      const Point_ref rpt = get(vpmap_, source(prev_he, mesh_));

      old_min_angle = (std::min)(old_min_angle,
                                 (std::min)(get_angle(Vector(old_pos, lpt), Vector(old_pos, rpt)),
                                            (std::min)(get_angle(Vector(lpt, rpt), Vector(lpt, old_pos)),
                                                       get_angle(Vector(rpt, old_pos), Vector(rpt, lpt)))));
    }
    std::cout << "old min angle: " << old_min_angle << std::endl;

    for(halfedge_descriptor main_he : halfedges_around_source(v, mesh_))
    {
      const halfedge_descriptor prev_he = prev(main_he, mesh_);
      const Point_ref lpt = get(vpmap_, target(main_he, mesh_));
      const Point_ref rpt = get(vpmap_, source(prev_he, mesh_));

      std::cout << "new angles: " << std::endl;
      std::cout << get_angle(Vector(new_pos, lpt), Vector(new_pos, rpt)) << " ";
      std::cout << get_angle(Vector(lpt, rpt), Vector(lpt, new_pos)) << " ";
      std::cout << get_angle(Vector(rpt, new_pos), Vector(rpt, lpt)) << std::endl;

      if(get_angle(Vector(new_pos, lpt), Vector(new_pos, rpt)) < old_min_angle)
        return false;
      if(get_angle(Vector(lpt, rpt), Vector(lpt, new_pos)) < old_min_angle)
        return false;
      if(get_angle(Vector(rpt, new_pos), Vector(rpt, lpt)) < old_min_angle)
        return false;
    }

    return true;
  }

  // gradient descent
  // ----------------
  bool gradient_descent(const vertex_descriptor v,
                        const double precision)
  {
    bool move_flag;
    double x, y, z, x_new, y_new, z_new, drdx, drdy, drdz;
    x = get(vpmap_, v).x();
    y = get(vpmap_, v).y();
    z = get(vpmap_, v).z();

    double S_av = compute_average_area_around(v);
    double energy = measure_energy(v, S_av);

    // if the adjacent areas are absolutely equal
    if(energy == 0.)
      return false;

    double energy_new = 0;
    double relative_energy = 1;
    unsigned int t = 1;
    double eta0 = 0.01;
    //double power_t = 0.25;
    double t0 = 0.001;
    double eta = eta0 / (1 + t0*t);

    while(relative_energy > precision)
    {
      drdx = 0.;
      drdy = 0.;
      drdz = 0.;
      compute_derivatives(drdx, drdy, drdz, v, S_av);

      x_new = x - eta * drdx;
      y_new = y - eta * drdy;
      z_new = z - eta * drdz;

      Point_3 moved(x_new, y_new, z_new);
      energy_new = measure_energy(v, S_av, moved);

      if(energy_new < energy)
      {
        put(vpmap_, v, moved);
        move_flag = true;
      }
      else
        return false;

      relative_energy = CGAL::to_double((energy - energy_new) / energy);

      // update
      x = x_new;
      y = y_new;
      z = z_new;
      energy = energy_new;
      ++t;

      // could use eta = eta0 / pow(t, power_t);
      eta = eta0 / (1 + t0 * t);
    }

    return move_flag;
  }

  double element_area(const vertex_descriptor v1,
                      const vertex_descriptor v2,
                      const vertex_descriptor v3) const
  {
    return CGAL::to_double(CGAL::approximate_sqrt(
                       traits_.compute_squared_area_3_object()(get(vpmap_, v1),
                                                               get(vpmap_, v2),
                                                               get(vpmap_, v3))));
  }

  double element_area(const Point_3& P,
                      const vertex_descriptor v2,
                      const vertex_descriptor v3) const
  {
    return CGAL::to_double(CGAL::approximate_sqrt(
                       traits_.compute_squared_area_3_object()(P,
                                                               get(vpmap_, v2),
                                                               get(vpmap_, v3))));
  }

  void compute_derivatives(double& drdx, double& drdy, double& drdz,
                           const vertex_descriptor v,
                           const double S_av)
  {
    for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
    {
      vertex_descriptor pi = source(next(h, mesh_), mesh_);
      vertex_descriptor pi1 = target(next(h, mesh_), mesh_);
      double S = element_area(v, pi, pi1);

      Vector vec(get(vpmap_, pi), get(vpmap_, pi1));

      // minimize r:
      // r = Σ(S-S_av)^2
      // dr/dx = 2 Σ(S - S_av) dS/dx
      // area of triangle with respect to (x_a, y_a, z_a) =
      // (1/2) [(v_z - v_y)x_a + (v_x - v_z)y_a + (v_y - v_x)z_a + constants]
      // vector v is (x_c - x_b, y_c - y_b, z_c - z_b)
      drdx += (S - S_av) * 0.5 * (vec.z() - vec.y());
      drdy += (S - S_av) * 0.5 * (vec.x() - vec.z());
      drdz += (S - S_av) * 0.5 * (vec.y() - vec.x());
    }

    drdx *= 2;
    drdy *= 2;
    drdz *= 2;
  }

  double compute_average_area_around(const vertex_descriptor v)
  {
    double sum_areas = 0.;
    unsigned int number_of_edges = 0;

    for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
    {
      // opposite vertices
      vertex_descriptor pi = source(next(h, mesh_), mesh_);
      vertex_descriptor pi1 = target(next(h, mesh_), mesh_);

      double S = element_area(v, pi, pi1);
      sum_areas += S;
      ++number_of_edges;
    }

    return sum_areas / number_of_edges;
  }

  double measure_energy(const vertex_descriptor v,
                        const double S_av)
  {
    double energy = 0;
    unsigned int number_of_edges = 0;

    for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
    {
      vertex_descriptor pi = source(next(h, mesh_), mesh_);
      vertex_descriptor pi1 = target(next(h, mesh_), mesh_);
      double S = element_area(v, pi, pi1);

      energy += (S - S_av)*(S - S_av);
      ++number_of_edges;
    }

    return to_double(energy / number_of_edges);
  }

  double measure_energy(const vertex_descriptor v,
                        const double S_av,
                        const Point_3& new_P)
  {
    double energy = 0;
    unsigned int number_of_edges = 0;
    for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
    {
      vertex_descriptor pi = source(next(h, mesh_), mesh_);
      vertex_descriptor pi1 = target(next(h, mesh_), mesh_);
      double S = element_area(new_P, pi, pi1);

      energy += (S - S_av)*(S - S_av);
      ++number_of_edges;
    }

    return to_double(energy / (2 * number_of_edges));
  }

  bool is_constrained(const vertex_descriptor v)
  {
    return get(vcmap_, v);
  }

  template<typename FaceRange>
  void set_vertex_range(const FaceRange& face_range)
  {
    // reserve 3 * #faces space
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

private:

  // data members
  // ------------
  PolygonMesh& mesh_;
  VertexPointMap& vpmap_;
  VertexConstraintMap vcmap_;
  GeomTraits traits_;

  Tree* tree_ptr_;
  std::vector<vertex_descriptor> vrange_;
};

} // namespace internal
} // namespace Polygon_mesh_processing
} // namespace CGAL


#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_MESH_SMOOTHING_IMPL_H
