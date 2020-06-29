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

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CURVATURE_FLOW_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CURVATURE_FLOW_IMPL_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/Weights.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/Dynamic_property_map.h>
#include <CGAL/utility.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>
#include <unordered_map>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

// Empirically, _Meyer seems to produce the best results from the various weights available in Weights.h
template<typename TriangleMesh,
         typename VertexPointMap,
         typename CotangentValue = CGAL::internal::Cotangent_value_Meyer<TriangleMesh, VertexPointMap> >
struct Edge_cotangent_weight
  : public CotangentValue
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;

  Edge_cotangent_weight(TriangleMesh& pmesh_, VertexPointMap vpmap_) : CotangentValue(pmesh_, vpmap_) {}

  TriangleMesh& pmesh() { return CotangentValue::pmesh(); }

  double operator()(halfedge_descriptor he)
  {
    if(is_border_edge(he, pmesh()))
    {
      halfedge_descriptor h1 = next(he, pmesh());
      vertex_descriptor vs = source(he, pmesh());
      vertex_descriptor vt = target(he, pmesh());
      vertex_descriptor v1 = target(h1, pmesh());

      return CotangentValue::operator()(vs, v1, vt);
    }
    else
    {
      halfedge_descriptor h1 = next(he, pmesh());
      halfedge_descriptor h2 = prev(opposite(he, pmesh()), pmesh());
      vertex_descriptor vs = source(he, pmesh());
      vertex_descriptor vt = target(he, pmesh());
      vertex_descriptor v1 = target(h1, pmesh());
      vertex_descriptor v2 = source(h2, pmesh());

      return CotangentValue::operator()(vs, v1, vt) + CotangentValue::operator()(vs, v2, vt);
    }
  }
};

template<typename TriangleMesh,
         typename VertexPointMap,
         typename VertexConstraintMap,
         typename SparseLinearSolver,
         typename GeomTraits>
class Shape_smoother
{
  typedef typename GeomTraits::FT                                                 FT;
  typedef typename GeomTraits::Point_3                                            Point;
  typedef typename GeomTraits::Vector_3                                           Vector;
  typedef typename boost::property_traits<VertexPointMap>::reference              Point_ref;

  typedef CGAL::Triple<std::size_t, std::size_t, double>                          Triplet;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor             edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor             face_descriptor;

  typedef CGAL::dynamic_vertex_property_t<std::size_t>                            Vertex_local_index;
  typedef typename boost::property_map<TriangleMesh, Vertex_local_index>::type    IndexMap;

  // linear system
  typedef typename SparseLinearSolver::Matrix                                     Eigen_matrix;
  typedef typename SparseLinearSolver::Vector                                     Eigen_vector;

public:
  Shape_smoother(TriangleMesh& mesh,
                 VertexPointMap& vpmap,
                 VertexConstraintMap& vcmap,
                 const GeomTraits& traits)
    :
      mesh_(mesh),
      vpmap_(vpmap),
      vcmap_(vcmap),
      vimap_(get(Vertex_local_index(), mesh_)),
      scale_volume_after_smoothing(true),
      traits_(traits),
      weight_calculator_(mesh, vpmap)
  { }

  template<typename FaceRange>
  void init_smoothing(const FaceRange& face_range)
  {
    set_face_range(face_range);

    std::size_t id = 0;
    for(vertex_descriptor v : vertices(mesh_))
      put(vimap_, v, id++);

    // vertices that are not in the range or are constrained still need value '1' in D because the RHS is D * X^n
    diagonal_.assign(vertices(mesh_).size(), 1.);
    constrained_flags_.assign(vertices(mesh_).size(), false);

    for(vertex_descriptor v : vertices(mesh_))
    {
      if(is_constrained(v))
      {
        constrained_flags_[get(vimap_, v)] = true;

        // scaling things cannot preserve the position of more than a single constrained point
        if(anchor_point == boost::none)
          anchor_point = get(vpmap_, v);
        else
          scale_volume_after_smoothing = false;
      }
    }

    if(!CGAL::is_closed(mesh_))
      scale_volume_after_smoothing = false;
  }

  void setup_system(Eigen_matrix& A,
                    Eigen_vector& bx, Eigen_vector& by, Eigen_vector& bz,
                    std::vector<Triplet>& stiffness_elements,
                    const double& time)
  {
    compute_coefficient_matrix(A, stiffness_elements, time);
    compute_rhs(bx, by, bz);
  }

  bool solve_system(const Eigen_matrix& A,
                    Eigen_vector& Xx, Eigen_vector& Xy, Eigen_vector& Xz,
                    const Eigen_vector& bx, const Eigen_vector& by, const Eigen_vector& bz,
                    SparseLinearSolver& solver)
  {
    FT D;

    // calls compute once to factorize with the preconditioner
    if(!solver.factor(A, D))
    {
#ifdef CGAL_PMP_SMOOTHING_DEBUG
      std::cerr << "Could not factorize linear system with preconditioner." << std::endl;
#endif
      return false;
    }

    if(!solver.linear_solver(bx, Xx) ||
       !solver.linear_solver(by, Xy) ||
       !solver.linear_solver(bz, Xz))
    {
#ifdef CGAL_PMP_SMOOTHING_DEBUG
      std::cerr << "Could not solve linear system." << std::endl;
#endif
      return false;
    }

    return true;
  }

  void calculate_stiffness_matrix_elements(std::vector<Triplet>& stiffness_elements)
  {
    CGAL_assertion(stiffness_elements.empty());
    stiffness_elements.reserve(8 * vrange_.size());

    std::unordered_map<std::size_t, double> diag_coeff;
    for(face_descriptor f : frange_)
    {
      for(halfedge_descriptor hi : halfedges_around_face(halfedge(f, mesh_), mesh_))
      {
        // Get a single canonical non-border halfedge per edge
        if(is_border(hi, mesh_))
          continue;

        const halfedge_descriptor hi_opp = opposite(hi, mesh_);
        if(!is_border(hi_opp, mesh_) && hi < hi_opp)
          continue;

        const vertex_descriptor v_source = source(hi, mesh_);
        const vertex_descriptor v_target = target(hi, mesh_);

        const bool is_source_constrained = is_constrained(v_source);
        const bool is_target_constrained = is_constrained(v_target);

        if(is_source_constrained && is_target_constrained)
          continue;

        const FT Lij = weight_calculator_(hi);

        const std::size_t i_source = get(vimap_, v_source);
        const std::size_t i_target = get(vimap_, v_target);

        // note that these constraints create asymmetry in the matrix
        if(!is_source_constrained)
        {
          stiffness_elements.push_back(Triplet(i_source, i_target, Lij));
          diag_coeff.insert(std::make_pair(i_source, 0)).first->second -= Lij;
        }

        if(!is_target_constrained)
        {
          stiffness_elements.push_back(Triplet(i_target, i_source, Lij));
          diag_coeff.insert(std::make_pair(i_target, 0)).first->second -= Lij;
        }
      }
    }

    typename std::unordered_map<std::size_t, double>::iterator it = diag_coeff.begin(),
                                                               end = diag_coeff.end();
    for(; it!=end; ++it)
      stiffness_elements.push_back(Triplet(it->first, it->first, it->second));
  }

  void update_mesh_no_scaling(const Eigen_vector& Xx, const Eigen_vector& Xy, const Eigen_vector& Xz)
  {
    for(vertex_descriptor v : vrange_)
    {
      std::size_t index = get(vimap_, v);
      const FT x_new = Xx[index];
      const FT y_new = Xy[index];
      const FT z_new = Xz[index];

      Point new_pos(x_new, y_new, z_new);
      put(vpmap_, v, new_pos);
    }
  }

  void update_mesh(const Eigen_vector& Xx, const Eigen_vector& Xy, const Eigen_vector& Xz)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    if(!scale_volume_after_smoothing)
      return update_mesh_no_scaling(Xx, Xy, Xz);

    const FT old_vol = volume(mesh_, parameters::vertex_point_map(vpmap_).geom_traits(traits_));

    // If no vertex is constrained, then the smoothed mesh will simply share the same centroid as the input mesh
    Point pre_smooth_anchor_point;
    if(anchor_point != boost::none)
      pre_smooth_anchor_point = *anchor_point;
    else
      pre_smooth_anchor_point = PMP::centroid(mesh_, parameters::vertex_point_map(vpmap_).geom_traits(traits_));

    for(vertex_descriptor v : vrange_)
    {
      std::size_t index = get(vimap_, v);
      const FT x_new = Xx[index];
      const FT y_new = Xy[index];
      const FT z_new = Xz[index];

      Point new_pos(x_new, y_new, z_new);
      put(vpmap_, v, new_pos);
    }

    Point post_smooth_anchor_point;
    if(anchor_point != boost::none)
      post_smooth_anchor_point = *anchor_point;
    else
      post_smooth_anchor_point = PMP::centroid(mesh_, parameters::vertex_point_map(vpmap_).geom_traits(traits_));

    const FT new_vol = volume(mesh_, parameters::vertex_point_map(vpmap_));
    CGAL_assertion(new_vol != 0);
    const FT inflating_factor = std::cbrt(CGAL::abs(old_vol / new_vol));

    for(vertex_descriptor v : vertices(mesh_))
    {
      Vector d = traits_.construct_vector_3_object()(post_smooth_anchor_point, get(vpmap_, v));
      Point new_pos = traits_.construct_translated_point_3_object()(pre_smooth_anchor_point,
                        traits_.construct_scaled_vector_3_object()(d, inflating_factor));

      put(vpmap_, v, new_pos);
    }
  }

private:
  bool is_constrained(const vertex_descriptor& v)
  {
    return get(vcmap_, v);
  }

  template<typename FaceRange>
  void set_face_range(const FaceRange& face_range)
  {
    frange_.assign(face_range.begin(), face_range.end());
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

  void compute_coefficient_matrix(Eigen_matrix& A,
                                  std::vector<Triplet>& stiffness_elements,
                                  const double& time)
  {
    fill_mass_matrix();

    // fill A = Mass - time * Laplacian
    for(const Triplet& t : stiffness_elements)
    {
      std::size_t i = t.get<0>(), j = t.get<1>();

      if(i != j)
        A.set_coef(i, j, - time * t.get<2>(), true);
      else if(!constrained_flags_[i]) // && i==j
        A.set_coef(i, i, diagonal_[t.get<0>()] - time * t.get<2>(), true);
    }

    for(vertex_descriptor v : vrange_)
    {
      std::size_t index = get(vimap_, v);
      if(constrained_flags_[index])
        A.set_coef(index, index, 1., true);
    }

    // we do not call A.assemble_matrix here
    // Eigen's compute during factorization does the building correctly,
    // and without assemble_matrix the reference A can be used in the next iterations.
  }

  void fill_mass_matrix()
  {
    for(vertex_descriptor v : vrange_)
    {
      std::size_t index = get(vimap_, v);
      if(!is_constrained(v))
        diagonal_[index] = 0.;
    }

    for(face_descriptor f : frange_)
    {
      const double area = face_area(f, mesh_, parameters::vertex_point_map(vpmap_).geom_traits(traits_));

      for(vertex_descriptor v : vertices_around_face(halfedge(f, mesh_), mesh_))
      {
        if(!is_constrained(v))
          diagonal_[get(vimap_, v)] += area / 6.;
      }
    }
  }

  void compute_rhs(Eigen_vector& bx, Eigen_vector& by, Eigen_vector& bz)
  {
    for(vertex_descriptor vi : vertices(mesh_))
    {
      std::size_t index = get(vimap_, vi);
      Point_ref p = get(vpmap_, vi);
      bx.set(index, diagonal_[index] * p.x());
      by.set(index, diagonal_[index] * p.y());
      bz.set(index, diagonal_[index] * p.z());
    }
  }

private:
  std::vector<vertex_descriptor> vrange_;
  std::vector<face_descriptor> frange_;
  TriangleMesh& mesh_;

  VertexPointMap vpmap_;
  VertexConstraintMap vcmap_;
  IndexMap vimap_;

  // Smoothing has a tendency to reduce volumes, so we can scale things back up based on the change
  // of volume. We need an anchor point to scale up, either a constrained point or the centroid
  // of the initial mesh if no vertex is constrained. If there is more than a constrained vertex,
  // then no scaling can be done without violating the constraint.
  bool scale_volume_after_smoothing;
  boost::optional<Point> anchor_point;

  // linear system data
  std::vector<double> diagonal_; // index of vector -> index of vimap_
  std::vector<bool> constrained_flags_;

  const GeomTraits& traits_;
  Edge_cotangent_weight<TriangleMesh, VertexPointMap> weight_calculator_;
};

} // internal
} // PMP
} // CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CURVATURE_FLOW_NEW_IMPL_H
