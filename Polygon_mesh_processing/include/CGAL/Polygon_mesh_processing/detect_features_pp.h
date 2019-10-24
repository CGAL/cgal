// Copyright (c) 2019 GeometryFactory (France).
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
//
#ifndef CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYGON_MESH_PP_H
#define CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYGON_MESH_PP_H

#include <CGAL/license/Polygon_mesh_processing/detect_features.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/utils.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_diagonalize_traits.h>
#endif
#include <CGAL/Diagonalize_traits.h>

#include <iostream>
#include <limits>
#include <map>

#include <sstream>
#include <fstream>

// @fixme compare scalar products instead of computing so many angles

// @fixme anchor_dense --> filter small salient curves with no strong DA ?

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

// Based on the paper
// Identification of C^1 and C^2 Discontinuities for Surface Meshes in CAD
// Xiangmin Jiao and Narasimha R. Bayyana

template <typename PolygonMesh, typename VPM, typename GeomTraits, typename EIFMap>
class Detector
{
  typedef typename GeomTraits::FT                                                 FT;
  typedef typename GeomTraits::Vector_3                                           Vector;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor            vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor          halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor              edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor              face_descriptor;

  typedef typename boost::property_traits<VPM>::value_type                        Point;
  typedef typename boost::property_traits<VPM>::reference                         Point_reference;

  typedef CGAL::dynamic_vertex_property_t<bool>                                   Vertex_bool_tag;
  typedef CGAL::dynamic_vertex_property_t<double>                                 Vertex_double_tag;
  typedef CGAL::dynamic_vertex_property_t<int>                                    Vertex_int_tag;
  typedef CGAL::dynamic_halfedge_property_t<double>                               Halfedge_double_tag;
  typedef CGAL::dynamic_edge_property_t<double>                                   Edge_double_tag;
  typedef CGAL::dynamic_edge_property_t<bool>                                     Edge_bool_tag;

  typedef typename boost::property_map<PolygonMesh, Vertex_bool_tag>::type        Vertex_bool_pmap;
  typedef typename boost::property_map<PolygonMesh, Vertex_double_tag>::type      Vertex_double_pmap;
  typedef typename boost::property_map<PolygonMesh, Vertex_int_tag>::type         Vertex_int_pmap;
  typedef typename boost::property_map<PolygonMesh, Halfedge_double_tag>::type    Halfedge_double_pmap;
  typedef typename boost::property_map<PolygonMesh, Edge_bool_tag>::type          Edge_bool_pmap;
  typedef typename boost::property_map<PolygonMesh, Edge_double_tag>::type        Edge_double_pmap;

  typedef CGAL::Eigen_diagonalize_traits<double, 3>                               Diagonalize_traits;
  typedef Eigen::Matrix3d                                                         Eigen_matrix;
  typedef Eigen::Vector3d                                                         Eigen_vector;

public:
  Detector(PolygonMesh& mesh,
           const VPM vpmap,
           const GeomTraits& traits,
           const EIFMap ifemap /*is_feature_edge_map*/)
    : mesh_(mesh), vpmap_(vpmap), gt_(traits), ifemap_(ifemap)
  { }

private:
  void print_candidate_edges(const std::string filename) const
  {
    std::ofstream out(filename);
    out.precision(17);

    for(edge_descriptor e : edges(mesh_))
    {
      if(get(candidate_edges_, e))
        out << "2 " << get(vpmap_, source(e, mesh_)) << " " << get(vpmap_, target(e, mesh_)) << std::endl;
    }
  }

  // -----------------------------------------------------------------------------------------------
  // ------------------------------------- ANGLE MEASURES ------------------------------------------
  // -----------------------------------------------------------------------------------------------
  double angle(const Point& p, const Point& q, const Point& r) const
  {
    return gt_.compute_approximate_angle_3_object()(p, q, r);
  }

  double angle_around_vertex(const vertex_descriptor v) const
  {
    double angle_sum = 0;
    for(halfedge_descriptor h : CGAL::halfedges_around_target(v, mesh_))
    {
      if(is_border(h, mesh_))
        continue;

      angle_sum += gt_.compute_approximate_angle_3_object()(get(vpmap_, source(h, mesh_)),
                                                            get(vpmap_, target(h, mesh_)),
                                                            get(vpmap_, target(next(h, mesh_), mesh_)));
    }

    return angle_sum;
  }

  double angle_defect(const vertex_descriptor v) const
  {
    const double total_angle = angle_around_vertex(v);
    if(is_border(v, mesh_))
      return 2 * (90 - total_angle);
    else
      return 2 * 90 - total_angle;
  }

  double dihedral_angle(const halfedge_descriptor h) const
  {
    if(is_border_edge(h, mesh_))
      return 90;

    Point_reference p = get(vpmap_, source(h, mesh_));
    Point_reference q = get(vpmap_, target(h, mesh_));
    Point_reference r = get(vpmap_, target(next(h, mesh_), mesh_));
    Point_reference s = get(vpmap_, target(next(opposite(h, mesh_), mesh_), mesh_));

    if(gt_.coplanar_3_object()(p, q, r, s))
      return 0.;

    // "Dihedral angle" is here the angle between the normals, which is supplementary of what is called
    // dihedral angle in the CGAL kernel...
    double val = gt_.compute_approximate_dihedral_angle_3_object()(p, q, r, s);

    val = (val >= 0) ? 180 - val : -180 - val;

    return val;
  }

  double dihedral_angle(const edge_descriptor e) const
  {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "dihedral_angle(" << e << ") -- [" << get(vpmap_, source(e, mesh_)) << "] [" << get(vpmap_, target(e, mesh_)) << "]" << std::endl;
    std::cout << "  dihedral_angle(" << halfedge(e, mesh_) << "): " << dihedral_angle(halfedge(e, mesh_)) << std::endl;
    std::cout << "  dihedral_angle(" << opposite(halfedge(e, mesh_), mesh_) << "): " << dihedral_angle(opposite(halfedge(e, mesh_), mesh_)) << std::endl;
#endif

    // can happen with one being 180 and the other -180
//    CGAL_warning(CGAL::abs(dihedral_angle(halfedge(e, mesh_)) -
//                           dihedral_angle(opposite(halfedge(e, mesh_), mesh_))) < 1e-12 );

    return dihedral_angle(halfedge(e, mesh_));
  }

  double turning_angle(const halfedge_descriptor first_h, const halfedge_descriptor second_h) const
  {
    return angle(get(vpmap_, target(first_h, mesh_)),
                 get(vpmap_, source(first_h, mesh_)),
                 get(vpmap_, target(second_h, mesh_)));
  }

  // obviously @cache that
  std::pair<Vector, bool> ridge_direction(const vertex_descriptor v) const
  {
    double eps = CGAL::square(std::tan(0.5 * theta_f_));

    // construct the matrix
    Eigen_matrix M = Eigen_matrix::Zero();

    for(halfedge_descriptor h : CGAL::halfedges_around_target(v, mesh_))
    {
      if(is_border(h, mesh_))
        continue;

      const face_descriptor f = face(h, mesh_);
      CGAL_assertion(f != boost::graph_traits<PolygonMesh>::null_face());

      const Vector fn = Polygon_mesh_processing::compute_face_normal(f, mesh_, CGAL::parameters::vertex_point_map(vpmap_)
                                                                                                .geom_traits(gt_));

      const double face_area = Polygon_mesh_processing::face_area(f, mesh_, CGAL::parameters::vertex_point_map(vpmap_)
                                                                                             .geom_traits(gt_));

      Eigen_vector mj;
      mj << fn[0], fn[1], fn[2];

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
//      std::cout << "mj: " << std::endl << mj << std::endl;
//      std::cout << "part matrix: " << std::endl << mj * mj.transpose() << std::endl;
#endif

      M += face_area * mj * mj.transpose();

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
//      std::cout << "cumulated: " << std::endl << M << std::endl;
#endif
    }

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
//    std::cout << "final: " << std::endl << M << std::endl;
#endif

    // compute the eigenvalues/vectors
    Eigen::SelfAdjointEigenSolver<Eigen_matrix> eigen_solver;

    eigen_solver.computeDirect(M);

    std::map<double, Vector> sorted_eigenvalues;
    for(int i=0; i<3; ++i)
    {
      sorted_eigenvalues.insert(std::make_pair(eigen_solver.eigenvalues()[i],
                                               gt_.construct_vector_3_object()(eigen_solver.eigenvectors()(0, i),
                                                                               eigen_solver.eigenvectors()(1, i),
                                                                               eigen_solver.eigenvectors()(2, i))));
    }

    typename std::map<double, Vector>::iterator it = sorted_eigenvalues.begin();

    bool is_usable = sorted_eigenvalues.size();
    if(is_usable)
    {
      const double lambda_1 = (it++)->first;
      const double lambda_2 = (it++)->first;
      const double lambda_3 = it->first;

      is_usable = (lambda_2 >= eps * lambda_1 && lambda_3 >= 0.7 * lambda_2);
    }

    if(is_usable)
      return std::make_pair(sorted_eigenvalues.begin()->second, true);
    else
      return std::make_pair(CGAL::NULL_VECTOR, false);
  }

  double one_side_turning_angle(const halfedge_descriptor h) const
  {
    const vertex_descriptor v = source(h, mesh_);
    const std::pair<Vector, bool> dir = ridge_direction(v);

    if(dir.second)
    {
      const Point r = gt_.construct_translated_point_3_object()(get(vpmap_, v), dir.first);
      return angle(get(vpmap_, target(h, mesh_)), get(vpmap_, v), r);
    }

    return 0.;
  }

  // -----------------------------------------------------------------------------------------------
  // ------------------------------------- OSTA/DA STRENGTH ----------------------------------------
  // -----------------------------------------------------------------------------------------------
  bool is_locally_strong_in_DA(const halfedge_descriptor h) const
  {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "is_locally_strong_in_DA: " << edge(h, mesh_) << " source(" << source(h, mesh_) << ")" << std::endl;
#endif

    if(is_border_edge(h, mesh_))
      return true;

    const edge_descriptor e = edge(h, mesh_);
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "dihedral angle: " << get(dihedral_angles_, e) << std::endl;
#endif

    // @todo might factorize all that jazz a bit better (a lot of 'get(dih_angles, e)')
    if(!is_f_strong_in_DA(e))
      return false;

    const double dih_angle = get(dihedral_angles_, e);
    const bool is_DA_positive = (dih_angle >= 0.);

    double max_other_DA_osta_smaller_than_half_pi = - std::numeric_limits<double>::max();
    double max_other_DA_osta_greater_than_half_pi = - std::numeric_limits<double>::max();

    for(halfedge_descriptor other_h : CGAL::halfedges_around_source(h, mesh_))
    {
      if(other_h == h)
        continue;

      const edge_descriptor other_e = edge(other_h, mesh_);

      if(!is_f_strong_in_DA(other_e))
        continue;

      const double other_dih_angle = get(dihedral_angles_, other_e);
      const bool is_other_DA_positive = (other_dih_angle >= 0.);

      if(is_DA_positive != is_other_DA_positive)
        continue;

      const double twice_other_osta = 2 * get(osta_angles_, other_h);

      if(twice_other_osta < 90) // <=> osta < pi/2
        max_other_DA_osta_smaller_than_half_pi = (CGAL::max)(max_other_DA_osta_smaller_than_half_pi,
                                                             other_dih_angle);
      else if(twice_other_osta >= 90) // <=> osta < pi/2
        max_other_DA_osta_greater_than_half_pi = (CGAL::max)(max_other_DA_osta_greater_than_half_pi,
                                                             other_dih_angle);
    }

    return (dih_angle >= max_other_DA_osta_smaller_than_half_pi ||
            dih_angle >= max_other_DA_osta_greater_than_half_pi);
  }

  bool is_locally_strong_in_OSTA(const halfedge_descriptor h) const
  {
    if(is_border_edge(h, mesh_))
      return true;

    const double osta_angle = get(osta_angles_, h);

    if(osta_angle <= theta_t_)
    {
      double min_osta_value = std::numeric_limits<double>::max();
      for(halfedge_descriptor other_h : CGAL::halfedges_around_source(h, mesh_)) // @fixme? border
      {
        if(other_h == h)
          continue;

        const edge_descriptor other_e = edge(other_h, mesh_);
        if(!is_f_strong_in_DA(other_e))
          continue;

        min_osta_value = (CGAL::min)(min_osta_value, get(osta_angles_, h));
      }

      if(osta_angle <= min_osta_value)
        return true;
    }
    else if(osta_angle >= 90 - theta_t_) // @strict
    {
      double max_osta_value = -(std::numeric_limits<double>::max());
      for(halfedge_descriptor other_h : CGAL::halfedges_around_source(h, mesh_)) // @fixme? border
      {
        if(other_h == h)
          continue;

        const edge_descriptor other_e = edge(other_h, mesh_);
        if(!is_f_strong_in_DA(other_e))
          continue;

        max_osta_value = (CGAL::max)(max_osta_value, get(osta_angles_, h));
      }

      if(osta_angle >= max_osta_value)
        return true;
    }

    return false;
  }

  bool is_e_strong_in_DA(const edge_descriptor e) const
  {
    return std::abs(get(dihedral_angles_, e)) >= theta_e_;
  }
  bool is_f_strong_in_DA(const edge_descriptor e) const
  {
    return std::abs(get(dihedral_angles_, e)) >= theta_f_;
  }
  bool is_k_strong_in_DA(const edge_descriptor e) const
  {
    return std::abs(get(dihedral_angles_, e)) >= theta_k_;
  }
  bool is_unconditionally_strong_in_DA(const edge_descriptor e) const
  {
    return std::abs(get(dihedral_angles_, e)) >= theta_F_;
  }
  bool is_unconditionally_strong_in_DA(const vertex_descriptor v) const
  {
    return get(angle_defects_, v) >= theta_D_;
  }

  bool is_sharp_corner(const vertex_descriptor v) const { return is_unconditionally_strong_in_DA(v); }
  bool is_sharp_edge(const edge_descriptor e) const { return is_unconditionally_strong_in_DA(e); }

  bool is_ambiguous_vertex(const vertex_descriptor v) const
  {
    // @cache
    std::pair<Vector, bool> dir = ridge_direction(v);

    return !dir.second;
  }

  // -----------------------------------------------------------------------------------------------
  // --------------------------------------- ATTACHEDNESS ------------------------------------------
  // -----------------------------------------------------------------------------------------------
  bool is_attached(const halfedge_descriptor h) const
  {
    // @fixme? weird formulation in the paper (p. 13 bottom) OR OR OR OR OR...

    const edge_descriptor e = edge(h, mesh_);
    if(!is_f_strong_in_DA(e))
      return false;

    const vertex_descriptor v = source(h, mesh_);

    for(halfedge_descriptor h : CGAL::halfedges_around_target(v, mesh_))
    {
      if(is_sharp_edge(edge(h, mesh_)))
        return true;
    }

    if(is_sharp_corner(v) || is_ambiguous_vertex(v))
      return true;

    if(is_locally_strong_in_OSTA(h) || is_locally_strong_in_DA(h))
      return true;

    return false;
  }

  bool is_strongly_attached(const halfedge_descriptor h) const
  {
    const edge_descriptor e = edge(h, mesh_);
    if(!is_f_strong_in_DA(e))
      return false;

    const vertex_descriptor v = source(h, mesh_);
    if(is_sharp_edge(edge(h, mesh_)) || is_sharp_corner(v) || is_ambiguous_vertex(v))
      return true;

    if(is_locally_strong_in_OSTA(h) && is_locally_strong_in_DA(h)) // note the '&&'
      return true;

    return false;
  }

  bool is_strongly_attached(const vertex_descriptor v) const
  {
    for(halfedge_descriptor h : CGAL::halfedges_around_source(v, mesh_))
    {
      if(is_strongly_attached(h))
        return true;
    }

    return false;
  }

  int count_incident_strongly_attached_halfedges(const vertex_descriptor v) const
  {
    int counter = 0;
    for(halfedge_descriptor h : CGAL::halfedges_around_source(v, mesh_))
    {
      if(is_strongly_attached(h))
        ++counter;
    }

    return counter;
  }

  bool is_quasi_strong(const halfedge_descriptor h) const
  {
    if(is_strongly_attached(source(h, mesh_)) &&
       is_strongly_attached(target(h, mesh_)) &&
       is_attached(h))
      return true;

    // protection for T-junctions
    if(count_incident_strongly_attached_halfedges(source(h, mesh_)) == 2 &&
       count_incident_strongly_attached_halfedges(target(h, mesh_)) == 2 &&
       is_locally_strong_in_OSTA(opposite(h, mesh_)) &&
       is_locally_strong_in_DA(opposite(h, mesh_)))
      return true;

    return false;
  }

  bool is_quasi_strong(const edge_descriptor e) const
  {
    halfedge_descriptor h = halfedge(e, mesh_);
    return (is_quasi_strong(h) && is_quasi_strong(opposite(h, mesh_)));
  }

  // -----------------------------------------------------------------------------------------------
  // ------------------------------------------- JOINTURE ------------------------------------------
  // -----------------------------------------------------------------------------------------------
  bool is_singleton(const halfedge_descriptor h) const
  {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "is " << edge(h, mesh_) << " source(" << source(h, mesh_) << ") singleton? ";
    std::cout << (get(candidate_edges_, edge(h, mesh_)) && get(candidate_valence_, source(h, mesh_)) == 1) << std::endl;
#endif

    return (get(candidate_edges_, edge(h, mesh_)) &&
            get(candidate_valence_, source(h, mesh_)) == 1);
  }

  bool is_dangling(const halfedge_descriptor h) const
  {
    if(!is_singleton(h))
      return false;

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "is " << edge(h, mesh_) << " source(" << source(h, mesh_) << ") dangling? ";
    std::cout << (!is_unconditionally_strong_in_DA(edge(h, mesh_)) || !is_sharp_corner(source(h, mesh_))) << std::endl;
#endif

    return (!is_unconditionally_strong_in_DA(edge(h, mesh_)) ||
            !is_sharp_corner(source(h, mesh_)));
  }

  bool is_semi_joint(const halfedge_descriptor h) const
  {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "is " << edge(h, mesh_) << " source(" << source(h, mesh_) << ") semi-joint? ";
#endif

    if(!get(candidate_edges_, edge(h, mesh_))) // @fixme ? not in paper, but can't see it being not requisite
      return false;

    // not using 'candidate_valence' because we care about the halfedges
    halfedge_descriptor first_h = h,
                        second_h = boost::graph_traits<PolygonMesh>::null_halfedge();
    for(halfedge_descriptor other_h : CGAL::halfedges_around_source(h, mesh_))
    {
      if(get(candidate_edges_, edge(other_h, mesh_)))
      {
        if(second_h == boost::graph_traits<PolygonMesh>::null_halfedge())
        {
          second_h = other_h;
        }
        else
        {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
          std::cout << "0" << std::endl;
#endif
          return false; // strictly more than 2 incident candidate halfedges
        }
      }
    }

    if(second_h == boost::graph_traits<PolygonMesh>::null_halfedge()) // only 1 incident candidate edge
    {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
      std::cout << "0" << std::endl;
#endif
      return false;
    }

    // Now we have exactly 2 candidate edges
    edge_descriptor first_e = edge(first_h, mesh_);
    edge_descriptor second_e = edge(second_h, mesh_);

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << ((is_sharp_corner(source(h, mesh_)) || turning_angle(first_h, second_h) > theta_T_) &&
                  (!is_unconditionally_strong_in_DA(first_e) || !is_unconditionally_strong_in_DA(second_e))) << std::endl;
#endif

    return ((is_sharp_corner(source(h, mesh_)) || turning_angle(first_h, second_h) > theta_T_) &&
            (!is_unconditionally_strong_in_DA(first_e) || !is_unconditionally_strong_in_DA(second_e)));
  }

  bool has_incident_unconditionally_strong_in_DA_incident_edges(const vertex_descriptor v) const
  {
    for(halfedge_descriptor h : CGAL::halfedges_around_target(v, mesh_))
    {
      if(is_unconditionally_strong_in_DA(edge(h, mesh_)))
        return true;
    }

    return false;
  }

  bool has_incident_acute_edge(const vertex_descriptor v) const
  {
    for(halfedge_descriptor h : CGAL::halfedges_around_target(v, mesh_))
    {
      if(is_border_edge(h, mesh_))
        continue;

      if(get(dihedral_angles_, edge(h, mesh_)) > 90) // @todo >= (check across the whole file actually)
        return true;
    }

    return false;
  }

  bool is_disjoint(const halfedge_descriptor h) const
  {
    const vertex_descriptor v = source(h, mesh_);

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "is " << edge(h, mesh_) << " source(" << source(h, mesh_) << ") disjoint?" << std::endl;
    std::cout << "valence: " << get(candidate_valence_, v) << std::endl;
#endif

    if(get(candidate_valence_, v) < 3)
      return false;

    const bool loc_strong_in_DA = is_locally_strong_in_DA(h);
    const bool loc_strong_in_OSTA = is_locally_strong_in_OSTA(h);
    const bool u_strong = is_unconditionally_strong_in_DA(edge(h, mesh_));
    const bool sharp = is_sharp_corner(v);
    const bool ambiguous = is_ambiguous_vertex(v);

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "is_locally_strong_in_DA: " << loc_strong_in_DA << std::endl;
    std::cout << "is_locally_strong_in_OSTA: " << loc_strong_in_OSTA << std::endl;
    std::cout << "is_unconditionally_strong_in_DA: " << u_strong << std::endl;
    std::cout << "is_sharp_corner: " << sharp << std::endl;
    std::cout << "is_ambiguous_vertex: " << ambiguous << std::endl;
#endif

    // @cache stuff (esp. 'is_locally_strong_in_...') and order from least expensive to most expensive
    const bool cond_1 = (!loc_strong_in_DA && !loc_strong_in_OSTA && !u_strong && !sharp && !ambiguous);
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "thus Disjoint cond_1: " << cond_1 << std::endl;
#endif
    if(cond_1)
      return true;

    const bool e_strong = is_e_strong_in_DA(edge(h, mesh_));
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "is_e_strong_in_DA: " << e_strong << std::endl;
#endif

    const bool cond_2 = ((!loc_strong_in_DA || !loc_strong_in_OSTA) && !e_strong && !sharp && !ambiguous);
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "thus Disjoint cond_2: " << cond_2 << std::endl;
#endif
    if(cond_2)
      return true;

    const bool has_incident_u_strong_edges = has_incident_unconditionally_strong_in_DA_incident_edges(v);

    const bool cond_3 = (!e_strong && has_incident_u_strong_edges);
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "thus Disjoint cond_3: " << cond_3 << std::endl;
#endif
    if(cond_3)
      return true;

    const bool cond_4 = (!u_strong && has_incident_acute_edge(v));
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "thus Disjoint cond_4: " << cond_4 << std::endl;
#endif
    if(cond_4)
      return true;

    return false;
  }

  bool is_multi_joint(const halfedge_descriptor h) const
  {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "is " << edge(h, mesh_) << " w/ source(" << source(h, mesh_) << ") multi-joint?" << std::endl;
    std::cout << "mj valence: " << get(candidate_valence_, source(h, mesh_)) << std::endl;
#endif

    return (get(candidate_valence_, source(h, mesh_)) >= 3 && !is_disjoint(h));
  }

  bool is_end_halfedge(const halfedge_descriptor h) const
  {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "is " << edge(h, mesh_) << " w/ source(" << source(h, mesh_) << ") an end?" << std::endl;
#endif

    // @cache all that stuff as much as possible
    return (is_singleton(h) || is_semi_joint(h) || is_disjoint(h) || is_multi_joint(h));
  }

  bool is_obscure_end_halfedge(const halfedge_descriptor h) const
  {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "is " << edge(h, mesh_) << " (source: " << source(h, mesh_) << " ) obscure? " << std::endl;
#endif

    return (is_dangling(h) || is_semi_joint(h) || is_disjoint(h));
  }

  bool is_salient_curve(const halfedge_descriptor start, const halfedge_descriptor last) const
  {
    if(target(last, mesh_) == source(start, mesh_)) // closed curve
      return true;

    if(is_obscure_end_halfedge(start) || is_obscure_end_halfedge(opposite(last, mesh_)))
      return false;

    // both end edges are non-obscure + not closed
    return true;
  }

  bool is_obscure_curve(const halfedge_descriptor start,
                        const halfedge_descriptor last,
                        const int number_of_edges_with_DAs_over_threshold) const
  {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "Is curve obscure? from " << edge(start, mesh_) << " to " << edge(last, mesh_) << std::endl;
#endif

    if(is_salient_curve(start, last))
      return false;

    const halfedge_descriptor last_opp = opposite(last, mesh_);
    const bool is_start_end_obscure = is_obscure_end_halfedge(start);
    const bool is_last_end_obscure = is_obscure_end_halfedge(last_opp);

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "start:end obscurity? " << is_start_end_obscure << " " << is_last_end_obscure << std::endl;
    std::cout << "start:end danglingness? " << is_dangling(start) << " " << is_dangling(last_opp) << std::endl;
    std::cout << "over DAs #: " << number_of_edges_with_DAs_over_threshold << std::endl;
#endif

    const bool cond_1 = ((is_start_end_obscure && is_last_end_obscure) ||
                         is_dangling(start) || is_dangling(last_opp)) && // @fixme xor?
                        number_of_edges_with_DAs_over_threshold < k_; // @fixme, with equality?
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "thus cond_1: " << cond_1 << std::endl;
#endif

    if(cond_1)
      return true;

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "unconditional DA (ex0): " << has_incident_unconditionally_strong_in_DA_incident_edges(source(start, mesh_)) << std::endl;
    std::cout << "unconditional DA (ex1): " << has_incident_unconditionally_strong_in_DA_incident_edges(target(last, mesh_)) << std::endl;
#endif

    const bool cond_2 = (is_start_end_obscure || is_last_end_obscure) && // @fixme xor?
                        has_incident_unconditionally_strong_in_DA_incident_edges(source(start, mesh_)) &&
                        has_incident_unconditionally_strong_in_DA_incident_edges(target(last, mesh_)) &&
                        number_of_edges_with_DAs_over_threshold == 0;
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "thus cond_2: " << cond_2 << std::endl;
#endif

    return cond_2;
  }

  bool is_semi_salient_curve(const halfedge_descriptor start, const halfedge_descriptor last) const
  {
    return !is_salient_curve(start, last) && !is_obscure_curve(start, last);
  }

  void tag_candidate_edges()
  {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "Tag candidate edges" << std::endl;
#endif

    for(vertex_descriptor v : vertices(mesh_))
      put(candidate_valence_, v, 0);

    int counter = 0; // @tmp
    for(edge_descriptor e : edges(mesh_))
    {
      const bool is_qs_edge = is_quasi_strong(e);
      put(candidate_edges_, e, is_qs_edge);

      if(is_qs_edge)
      {
        ++counter;
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
        std::cout << e << " -- [" << get(vpmap_, source(e, mesh_)) << "] ["
                                  << get(vpmap_, target(e, mesh_)) << "] is quasi strong" << std::endl;
#endif

        put(candidate_valence_, source(e, mesh_), get(candidate_valence_, source(e, mesh_)) + 1);
        put(candidate_valence_, target(e, mesh_), get(candidate_valence_, target(e, mesh_)) + 1);
      }
    }

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << counter << " candidate edges" << std::endl;
    print_candidate_edges("results/all_unfiltered_sharp_edges.polylines.txt");

    for(vertex_descriptor v : vertices(mesh_))
      std::cout << "valence[" << v << " -- (" << get(vpmap_, v) << ")] = " << get(candidate_valence_, v) << std::endl;
#endif
  }

  // @fixme completely ignoring the fact that sharp corners are not marked as candidates currently
  // (don't really care since we want to tag sharp edges, but might want to be able to tag
  // sharp vertices too).

  halfedge_descriptor next_candidate_edge(const halfedge_descriptor in_h) const
  {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "walking through: " << target(in_h, mesh_) << " valence " << get(candidate_valence_, target(in_h, mesh_)) << std::endl;
#endif
    CGAL_precondition(get(candidate_valence_, target(in_h, mesh_)) == 2);

    const edge_descriptor in_e = edge(in_h, mesh_);
    for(halfedge_descriptor h : CGAL::halfedges_around_source(target(in_h, mesh_), mesh_))
    {
      if(edge(h, mesh_) == in_e)
        continue;

      if(get(candidate_edges_, edge(h, mesh_)))
      {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
        std::cout << "onto " << edge(h, mesh_) << std::endl;
#endif
        return h;
      }
    }

    CGAL_assertion(false);
    return boost::graph_traits<PolygonMesh>::null_halfedge();
  }

  void walk_curve(const halfedge_descriptor start_h,
                  std::list<halfedge_descriptor>& curve_halfedges,
                  int& number_of_edges_with_DAs_over_threshold) const
  {
    halfedge_descriptor curr_h = start_h;
    for(;;)
    {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
      std::cout << "curr_h: " << edge(curr_h, mesh_) << " ["
                              << get(vpmap_, source(curr_h, mesh_)) << " "
                              << get(vpmap_, target(curr_h, mesh_)) << "]" << std::endl;
#endif

      CGAL_assertion(get(candidate_edges_, edge(curr_h, mesh_)));
      curve_halfedges.push_back(curr_h);

      edge_descriptor curr_e = edge(curr_h, mesh_);

      if(is_k_strong_in_DA(curr_e))
        ++number_of_edges_with_DAs_over_threshold;

      // is this correct?
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
      std::cout << "check next for curve end (" << edge(opposite(curr_h, mesh_), mesh_) << ")" << std::endl;
#endif
      if(is_end_halfedge(opposite(curr_h, mesh_)))
      {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
        std::cout << "opposite is the end of things" << std::endl;
#endif
        break;
      }

      CGAL_assertion_code(halfedge_descriptor old_h = curr_h;)
      curr_h = next_candidate_edge(curr_h);
      CGAL_assertion(old_h != curr_h);
      CGAL_assertion(target(old_h, mesh_) == source(curr_h, mesh_));

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
      std::cout << "check next for curve end (" << edge(curr_h, mesh_) << ")" << std::endl;
#endif
      if(is_end_halfedge(curr_h)) // meaning the next one is _not_ added to the walk
      {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
        std::cout << "next is the end of things" << std::endl;
#endif
        break;
      }

      // might run again into 'start_h', which might have been an end edge when end edges were collected,
      // but is not an end edge anymore (e.g. because an incident edge is not candidate anymore)
      if(curr_h == start_h)
        break;
    }
  }

  void initialize_bounds(const FT theta_F, const FT theta_f, const FT theta_D, const FT theta_T,
                         const FT theta_t, const FT theta_e, const FT theta_k, int k)
  {
    std::cout << "Initialize bounds" << std::endl;
    std::cout << "theta_F: " << theta_F << std::endl;
    std::cout << "theta_f: " << theta_f << std::endl;
    std::cout << "theta_D: " << theta_D << std::endl;
    std::cout << "theta_T: " << theta_T << std::endl;
    std::cout << "theta_t: " << theta_t << std::endl;
    std::cout << "theta_e: " << theta_e << std::endl;
    std::cout << "theta_k: " << theta_k << std::endl;
    std::cout << "k: " << k << std::endl;

    theta_F_ = theta_F;
    theta_f_ = theta_f;
    theta_D_ = theta_D;
    theta_T_ = theta_T;
    theta_t_ = theta_t;
    theta_e_ = theta_e;
    theta_k_ = theta_k;
    k_ = k;

//    CGAL_warning(theta_F_ >= 10 * theta_f_);
//    CGAL_warning(theta_T_ >= 2 * theta_t_);
//    CGAL_warning(theta_f_ <= theta_e_ && theta_e_ <= theta_F_);
  }

  void initialize_maps()
  {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "Initialize maps" << std::endl;
#endif

    angle_defects_ = get(Vertex_double_tag(), mesh_);
    osta_angles_ = get(Halfedge_double_tag(), mesh_);
    dihedral_angles_ = get(Edge_double_tag(), mesh_);

    candidate_valence_ = get(Vertex_int_tag(), mesh_);
    candidate_edges_ = get(Edge_bool_tag(), mesh_);

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "Compute angle defects..." << std::endl;
#endif
    for(vertex_descriptor v : vertices(mesh_))
      put(angle_defects_, v, angle_defect(v));

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "Compute dihedral angles..." << std::endl;
#endif
    for(edge_descriptor e : edges(mesh_))
      put(dihedral_angles_, e, dihedral_angle(e));

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "Compute OSTAs..." << std::endl;
#endif
    for(halfedge_descriptor h : halfedges(mesh_))
    {
      put(osta_angles_, h, one_side_turning_angle(h));

      // @cache local extremas at vertices (don't forget about the DA filter)
    }
  }

public:
  void tag_sharp_edges(const FT theta_F = 65, // all angles in degrees
                       const FT theta_f = 10,
                       const FT theta_D = 60,
                       const FT theta_T = 40,
                       const FT theta_t = 20,
                       const FT theta_e = 25,
                       const FT theta_k = 50,
                       int k = 5)
  {
    std::cout << "Tagging sharp edges..." << std::endl;

    initialize_bounds(theta_F, theta_f, theta_D, theta_T, theta_t, theta_e, theta_k, k);
    initialize_maps();

    // build list of incident candidate halfedges at each vertex
    tag_candidate_edges();

    int iter = 0;

    bool removed_some_halfedges;
    do
    {
      removed_some_halfedges = false;

      std::stringstream oss_i;
      oss_i << "results/iter_" << iter++ << "_sharp_edges.polylines.txt" << std::ends;
      print_candidate_edges(oss_i.str());

      // collect obscure end-edges from quasi-strong edges
      std::unordered_set<halfedge_descriptor> obscure_end_halfedges;
      for(halfedge_descriptor h : halfedges(mesh_))
      {
        // @cache
        if(get(candidate_edges_, edge(h, mesh_)) && is_obscure_end_halfedge(h))
          obscure_end_halfedges.insert(h);
      }

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
      std::cout << obscure_end_halfedges.size() << " obscure end halfedges" << std::endl;
#endif

      // while there are still obscure end-edges left
      while(!obscure_end_halfedges.empty())
      {
        // traverse curve
        halfedge_descriptor start_h = *(obscure_end_halfedges.begin());
        obscure_end_halfedges.erase(obscure_end_halfedges.begin());

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
        std::cout << "-----" << std::endl;
        std::cout << "start at: " << edge(start_h, mesh_) << " ["
                                  << get(vpmap_, source(start_h, mesh_)) << " "
                                  << get(vpmap_, target(start_h, mesh_)) << "]" << std::endl;
#endif

        // check if this is still a candidate or it's been cleaned up
        // simpler than removing the end edge of obscure curves (which wouldn't be enough in any case:
        // imagine a triple point with one edge being uncandidated, then surviving start edges
        // can now be in the middle of an obscure curve)
        if(!get(candidate_edges_, edge(start_h, mesh_)))
           continue;

        std::list<halfedge_descriptor> curve_halfedges;
        int number_of_edges_with_DAs_over_threshold = 0;

        walk_curve(start_h, curve_halfedges, number_of_edges_with_DAs_over_threshold);
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
        std::cout << "~~~~~> Curve has length: " << curve_halfedges.size() << std::endl;
#endif

        const halfedge_descriptor last_h = curve_halfedges.back();
        // check if the curve is obscure
        if(is_obscure_curve(start_h, last_h, number_of_edges_with_DAs_over_threshold))
        {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
          std::cout << "Obscure curve!!" << std::endl;
#endif

          for(halfedge_descriptor h : curve_halfedges)
          {
            put(candidate_edges_, edge(h, mesh_), false);
            CGAL_assertion(get(candidate_valence_, source(h, mesh_)) > 0);
            put(candidate_valence_, source(h, mesh_), get(candidate_valence_, source(h, mesh_)) - 1);
          }

          CGAL_assertion(get(candidate_valence_, target(last_h, mesh_)) > 0);
          put(candidate_valence_, target(last_h, mesh_), get(candidate_valence_, target(last_h, mesh_)) - 1);

          removed_some_halfedges = true;
        }
        else
        {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
          std::cout << "Curve is NOT obscure!!" << std::endl;
#endif
        }
      }
    }
    while(removed_some_halfedges); // while removed halfedges

#ifndef SHARP_EDGES_IN_INDEPENDENT_FILES
    std::ofstream out("results/all_sharp_edges.polylines.txt");
#endif

    // mark remaining candidates halfedges as C1 discontinuities
    int counter = 0; // @tmp
    for(edge_descriptor e : edges(mesh_))
    {
      const bool is_sharp = get(candidate_edges_, e);
      put(ifemap_, e, is_sharp);

      if(is_sharp)
      {
#ifdef SHARP_EDGES_IN_INDEPENDENT_FILES
        std::stringstream oss;
        oss << "results/sharp_edge_" << counter << ".polylines.txt" << std::ends;
        std::ofstream out(oss.str().c_str());
#endif
        out << "2 " << get(vpmap_, source(e, mesh_)) << " " << get(vpmap_, target(e, mesh_)) << std::endl;
        ++counter;
      }
    }

    std::cout << counter << " sharp edges" << std::endl;
  }

private:
  PolygonMesh& mesh_;
  const VPM vpmap_;
  GeomTraits gt_;
  const EIFMap ifemap_;

  Vertex_double_pmap angle_defects_;
  Halfedge_double_pmap osta_angles_;
  Edge_double_pmap dihedral_angles_;

  Vertex_int_pmap candidate_valence_;
  Edge_bool_pmap candidate_edges_;

  // bounds
  double theta_D_;
  double theta_e_;
  double theta_f_; // Bound to define strong edges (DA)
  double theta_F_; // Bound to define unconditionally strong edges (DA)
  double theta_k_;
  double theta_t_;
  double theta_T_;

  int k_;
};

} // end namespace internal

template <typename FaceRange,
          typename PolygonMesh,
          typename FT,
          typename EdgeIsFeatureMap,
          typename NamedParameters>
void detect_sharp_edges_pp(const FaceRange& /*faces*/,
                           PolygonMesh& pmesh,
                           const FT strong_DA_in_deg,
                           EdgeIsFeatureMap edge_is_feature_map,
                           const NamedParameters& np)
{
  using boost::choose_param;
  using boost::get_param;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type                VPM;
  VPM vpm = choose_param(get_param(np, internal_np::vertex_point),
                         get_const_property_map(CGAL::vertex_point, pmesh));

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type                    GeomTraits;
  GeomTraits traits = choose_param(get_param(np, internal_np::geom_traits), GeomTraits());

  typedef typename internal::Detector<PolygonMesh, VPM, GeomTraits, EdgeIsFeatureMap>   Detector;
  Detector detector(pmesh, vpm, traits, edge_is_feature_map);

  const double weak_DA_in_deg = choose_param(get_param(np, internal_np::weak_dihedral_angle), 10.);

  detector.tag_sharp_edges(strong_DA_in_deg, weak_DA_in_deg);
}

template <typename PolygonMesh,
          typename FT,
          typename EdgeIsFeatureMap,
          typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
void detect_sharp_edges_pp(PolygonMesh& pmesh,
                           const FT strong_DA_in_deg,
                           EdgeIsFeatureMap edge_is_feature_map,
                           const CGAL_PMP_NP_CLASS& np)
{
  return detect_sharp_edges_pp(faces(pmesh), pmesh, strong_DA_in_deg, edge_is_feature_map, np);
}

} // end namespace PMP
} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYGON_MESH_PP_H
