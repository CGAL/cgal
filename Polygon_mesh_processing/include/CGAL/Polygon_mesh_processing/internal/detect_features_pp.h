// Copyright (c) 2019 GeometryFactory (France).
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
//
#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_DETECT_FEATURES_IN_POLYGON_MESH_PP_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_DETECT_FEATURES_IN_POLYGON_MESH_PP_H

#include <CGAL/license/Polygon_mesh_processing/detect_features.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Named_function_parameters.h>
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

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

// Based on the paper
// Identification of C^1 and C^2 Discontinuities for Surface Meshes in CAD
// Xiangmin Jiao and Narasimha R. Bayyana

// The idea is to start from a large set of "candidate" edges and trim to a set of "sharp" edges
// using a large set of heurestics.

// Experimental code, probably has some bugs and improvements

// Usage: API gives access to two of the (many) parameters, so called "strong" and "weak" dihedral angles.
// Candidate edges have a dihedral angle greater than the weak dihedral angle.
// A dihedral angle greater than the strong value ensures that the edge will be marked as sharp.
// The same behavior as the basic CGAL::detect_features() can be obtained by choosing weak == strong

// @fixme compare scalar products instead of computing so many angles
// @fixme anchor_dense --> filter small salient curves with no strong DA ?

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
      return 2. * (90. - total_angle);
    else
      return 2. * 90. - total_angle;
  }

  double dihedral_angle(const halfedge_descriptor h) const
  {
    if(is_border_edge(h, mesh_))
      return 90.;

    Point_reference p = get(vpmap_, source(h, mesh_));
    Point_reference q = get(vpmap_, target(h, mesh_));
    Point_reference r = get(vpmap_, target(next(h, mesh_), mesh_));
    Point_reference s = get(vpmap_, target(next(opposite(h, mesh_), mesh_), mesh_));

    if(gt_.coplanar_3_object()(p, q, r, s))
      return 0.;

    // "Dihedral angle" is here the angle between the normals, which is the supplementary angle
    // of what is called computed in the CGAL functor (the angle between the planes)...
    double val = gt_.compute_approximate_dihedral_angle_3_object()(p, q, r, s);

    val = (val >= 0.) ? 180. - val : -180. - val;

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
    const double eps = CGAL::square(std::tan(0.5 * theta_f_));

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

    eigen_solver.compute(M);

    std::map<double, Vector> sorted_eigenvalues;
    for(int i=0; i<3; ++i)
    {
      CGAL_assertion(eigen_solver.eigenvalues()[i] >= - 100 * std::numeric_limits<double>::epsilon());
      sorted_eigenvalues.insert(std::make_pair(eigen_solver.eigenvalues()[i],
                                               gt_.construct_vector_3_object()(eigen_solver.eigenvectors()(0, i),
                                                                               eigen_solver.eigenvectors()(1, i),
                                                                               eigen_solver.eigenvectors()(2, i))));
    }

    typename std::map<double, Vector>::iterator it = sorted_eigenvalues.begin();

    bool is_usable = (sorted_eigenvalues.size() == 3);
    if(is_usable)
    {
      const double lambda_3 = (it++)->first;
      const double lambda_2 = (it++)->first;
      const double lambda_1 = it->first;

      is_usable = (lambda_2 >= eps * lambda_1 && lambda_3 >= 0.7 * lambda_2);
    }

    if(is_usable)
      return std::make_pair(sorted_eigenvalues.rbegin()->second, true);
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
    std::cout << "dihedral angle ["  << e << "] = " << get(dihedral_angles_, e) << std::endl;
#endif

    // @todo might factorize all that jazz a bit better (a lot of 'get(dih_angles, e)')
    if(!is_f_strong_in_DA(e))
      return false;

    const double dih_angle = get(dihedral_angles_, e);
    const bool is_DA_positive = (dih_angle >= 0.);

    double max_other_DA_osta_smaller_than_half_pi = std::numeric_limits<double>::lowest();
    double max_other_DA_osta_greater_than_half_pi = std::numeric_limits<double>::lowest();

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

      const double twice_other_osta = 2. * get(osta_angles_, other_h);

      if(twice_other_osta < 90.) // <=> osta < pi/2
        max_other_DA_osta_smaller_than_half_pi = (CGAL::max)(max_other_DA_osta_smaller_than_half_pi,
                                                             other_dih_angle);
      else if(twice_other_osta >= 90.) // <=> osta >= pi/2
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
      for(halfedge_descriptor other_h : CGAL::halfedges_around_source(h, mesh_))
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
    else if(osta_angle >= 90. - theta_t_) // @strict
    {
      double max_osta_value = std::numeric_limits<double>::lowest();
      for(halfedge_descriptor other_h : CGAL::halfedges_around_source(h, mesh_))
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
    const double dh = std::abs(get(dihedral_angles_, e));
    return (dh >= theta_e_ && (!use_upper_bound_on_DA_ || dh <= 180. - theta_e_));
  }
  bool is_f_strong_in_DA(const edge_descriptor e) const
  {
    const double dh = std::abs(get(dihedral_angles_, e));
    return (dh >= theta_f_ && (!use_upper_bound_on_DA_ || dh <= 180. - theta_f_));
  }
  bool is_k_strong_in_DA(const edge_descriptor e) const
  {
    const double dh = std::abs(get(dihedral_angles_, e));
    return (dh >= theta_k_ && (!use_upper_bound_on_DA_ || dh <= 180. - theta_k_));
  }
  bool is_unconditionally_strong_in_DA(const edge_descriptor e) const
  {
    const double dh = std::abs(get(dihedral_angles_, e));
    return (dh >= theta_F_ && (!use_upper_bound_on_DA_ || dh <= 180. - theta_F_));
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
      if(other_h != h && get(candidate_edges_, edge(other_h, mesh_)))
      {
        if(second_h == boost::graph_traits<PolygonMesh>::null_halfedge())
        {
          second_h = other_h;
        }
        else
        {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
          std::cout << "more than 2 incident candidate halfedges" << std::endl;
#endif
          return false; // strictly more than 2 incident candidate halfedges
        }
      }
    }

    if(second_h == boost::graph_traits<PolygonMesh>::null_halfedge()) // only 1 incident candidate edge
    {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
      std::cout << "only 1 incident candidate edge" << std::endl;
#endif
      return false;
    }

    // Now we have exactly 2 candidate edges
    edge_descriptor first_e = edge(first_h, mesh_);
    edge_descriptor second_e = edge(second_h, mesh_);

    const bool cond_1 = is_sharp_corner(source(h, mesh_)) || turning_angle(first_h, second_h) <= theta_T_;
    const bool cond_2 = !is_unconditionally_strong_in_DA(first_e) || !is_unconditionally_strong_in_DA(second_e);

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "valence 2, cond_1/2: " << cond_1 << " " << cond_2 << std::endl;
    std::cout << "sharp: " << is_sharp_corner(source(h, mesh_)) << std::endl;
    std::cout << "turning angle: " << turning_angle(first_h, second_h) << std::endl;
    std::cout << "u_DA_1: " << is_unconditionally_strong_in_DA(first_e)
              << " u_DA_2: " << is_unconditionally_strong_in_DA(second_e) << std::endl;
#endif

    return (cond_1 && cond_2);
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
    std::cout << "- is " << edge(h, mesh_) << " source(" << source(h, mesh_) << ") multi-joint?" << std::endl;
    std::cout << "mj valence: " << get(candidate_valence_, source(h, mesh_)) << std::endl;
#endif

    return (get(candidate_valence_, source(h, mesh_)) >= 3 && !is_disjoint(h));
  }

  bool is_end_halfedge(const halfedge_descriptor h) const
  {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "- is " << edge(h, mesh_) << " source(" << source(h, mesh_) << ") an end?" << std::endl;
#endif

    // @cache all that stuff as much as possible
    return (is_singleton(h) || is_semi_joint(h) || is_disjoint(h) || is_multi_joint(h));
  }

  bool is_obscure_end_halfedge(const halfedge_descriptor h) const
  {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "- is " << edge(h, mesh_) << " source(" << source(h, mesh_) << ") an obscure end?" << std::endl;
#endif

    return (is_dangling(h) || is_semi_joint(h) || is_disjoint(h));
  }

  bool is_salient_curve(const halfedge_descriptor start_h,
                        const halfedge_descriptor last_h) const
  {


    // both end edges are non-obscure
    return true;
  }

  bool is_negligible_curve(const std::list<halfedge_descriptor>& curve_halfedges,
                           const double feature_length) const
  {
//    const halfedge_descriptor start_h = curve_halfedges.front();
//    const halfedge_descriptor last_h = curve_halfedges.back();
//    bool is_start_h_incident = is_semi_joint(start_h) || is_multi_joint(start_h);
//    bool is_last_h_incident = is_semi_joint(last_h) || is_multi_joint(last_h);

    // @tentative Below is done not to remove small features incident to other features,
    // with the imaginary case of a slim paralleliped. Keep?
//    std::cout << "Incidence: " << is_start_h_incident << " " << is_last_h_incident << std::endl;
//    if(is_start_h_incident || is_last_h_incident)
//      continue;

    double curve_length = 0.;
    for(halfedge_descriptor h : curve_halfedges)
    {
      curve_length += CGAL::sqrt(CGAL::squared_distance(get(vpmap_, source(h, mesh_)),
                                                        get(vpmap_, target(h, mesh_))));
    }

    return curve_length <= feature_length;
  }

  bool is_obscure_curve(const std::list<halfedge_descriptor>& curve_halfedges,
                        const int number_of_edges_with_DAs_over_threshold) const
  {
    const halfedge_descriptor start_h = curve_halfedges.front();
    const halfedge_descriptor last_h = curve_halfedges.back();

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "Is curve obscure? from " << edge(start_h, mesh_) << " to " << edge(last_h, mesh_) << std::endl;
#endif

    const bool is_closed_curve = (target(last_h, mesh_) == source(start_h, mesh_));
    if(is_closed_curve)
      return false;

    const halfedge_descriptor last_h_opp = opposite(last_h, mesh_);
    const bool is_start_h_end_obscure = is_obscure_end_halfedge(start_h);
    const bool is_last_h_end_obscure = is_obscure_end_halfedge(last_h_opp);
    const bool is_start_h_dangling = is_dangling(start_h);
    const bool is_last_h_dangling = is_dangling(last_h_opp);

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "start/last: obscurity? " << is_start_h_end_obscure << " " << is_last_h_end_obscure << std::endl;
    std::cout << "start/last: danglingness? " << is_start_h_dangling << " " << is_last_h_dangling << std::endl;
    std::cout << "over DAs #: " << number_of_edges_with_DAs_over_threshold << std::endl;
#endif

    // If neither ends are obscure, the curve is salient and it is not obscure
    if(!is_start_h_end_obscure && !is_last_h_end_obscure)
      return false;

    const bool cond_1 = ((is_start_h_end_obscure && is_last_h_end_obscure) ||
                         is_start_h_dangling || is_last_h_dangling) && // @fixme xor?
                        number_of_edges_with_DAs_over_threshold <= k_;
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "thus cond_1: " << cond_1 << std::endl;
#endif

    if(cond_1)
      return true;

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "unconditional DA (ex0): " << has_incident_unconditionally_strong_in_DA_incident_edges(source(start_h, mesh_)) << std::endl;
    std::cout << "unconditional DA (ex1): " << has_incident_unconditionally_strong_in_DA_incident_edges(target(last_h, mesh_)) << std::endl;
#endif

    const bool cond_2 = (is_start_h_end_obscure || is_last_h_end_obscure) && // @fixme xor?
                        has_incident_unconditionally_strong_in_DA_incident_edges(source(start_h, mesh_)) &&
                        has_incident_unconditionally_strong_in_DA_incident_edges(target(last_h, mesh_)) &&
                        number_of_edges_with_DAs_over_threshold == 0;
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "thus cond_2: " << cond_2 << std::endl;
#endif

    return cond_2;
  }

  bool is_semi_salient_curve(const halfedge_descriptor start_h,
                             const halfedge_descriptor last_h) const
  {
    // @cache
    return !is_salient_curve(start_h, last_h) && !is_obscure_curve(start_h, last_h);
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
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
      std::cout << std::endl << "---- Evaluate edge " << e << std::endl;
#endif

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
    std::cout << counter << " initial candidate edges" << std::endl;
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

  halfedge_descriptor prev_candidate_edge(const halfedge_descriptor in_h) const
  {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << "walking (backwwards) through: " << source(in_h, mesh_) << " valence " << get(candidate_valence_, source(in_h, mesh_)) << std::endl;
#endif
    CGAL_precondition(get(candidate_valence_, source(in_h, mesh_)) == 2);

    const edge_descriptor in_e = edge(in_h, mesh_);
    for(halfedge_descriptor h : CGAL::halfedges_around_target(source(in_h, mesh_), mesh_))
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

  void untag_curve(const std::list<halfedge_descriptor>& curve_halfedges)
  {
    for(halfedge_descriptor h : curve_halfedges)
    {
      put(candidate_edges_, edge(h, mesh_), false);
      CGAL_assertion(get(candidate_valence_, source(h, mesh_)) > 0);
      put(candidate_valence_, source(h, mesh_), get(candidate_valence_, source(h, mesh_)) - 1);
    }

    const halfedge_descriptor last_h = curve_halfedges.back();

    CGAL_assertion(get(candidate_valence_, target(last_h, mesh_)) > 0);
    put(candidate_valence_, target(last_h, mesh_), get(candidate_valence_, target(last_h, mesh_)) - 1);
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

  void trim_candidate_list(const double feature_length)
  {
    int iter = 0;

    bool removed_some_halfedges;
    do
    {
      removed_some_halfedges = false;

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
    std::cout << std::endl << "Iteration: " << iter << std::endl;
#endif

      std::stringstream oss_i;
      oss_i << "results/iter_" << iter++ << "_sharp_edges.polylines.txt" << std::ends;
      print_candidate_edges(oss_i.str());

      // collect obscure end-edges from quasi-strong edges
      std::unordered_set<halfedge_descriptor> obscure_end_halfedges;
      for(halfedge_descriptor h : halfedges(mesh_))
      {
        // @cache
        // @tentative check if the halfedge belongs to a small cycle (e.g. size <= k_ ?), and if so, untag the full cycle
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

        // check if the curve is obscure
        if(is_negligible_curve(curve_halfedges, feature_length) ||
           is_obscure_curve(curve_halfedges, number_of_edges_with_DAs_over_threshold))
        {
#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
          std::cout << "Obscure curve!!" << std::endl;
#endif

          untag_curve(curve_halfedges);

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

    print_candidate_edges("results/sharp_edges_after_trimming.polylines.txt");
  }

  // Curves with length <= feature_length are removed
  //
  // @tentative, also bound the number of edges? Immediately lose some edges on e.g. fandisk that are
  // curves of length '3' and actual sharp edges
  void remove_small_features(const double feature_length)
  {
    std::cout << "Remove small cycles" << std::endl;

    Edge_bool_pmap visited_edges = get(Edge_bool_tag(), mesh_);

    // starts by tagging all curves with ends (non cycles)
    for(halfedge_descriptor start_h : halfedges(mesh_))
    {
      if(!get(candidate_edges_, edge(start_h, mesh_)) ||
         get(visited_edges, edge(start_h, mesh_)) ||
         !is_end_halfedge(start_h))
        continue;

      std::cout << "Start a new curve at " << edge(start_h, mesh_) << std::endl;

      std::list<halfedge_descriptor> curve_halfedges;
      int number_of_edges_with_DAs_over_threshold = 0;
      walk_curve(start_h, curve_halfedges, number_of_edges_with_DAs_over_threshold);

      double curve_length = 0.;
      for(const halfedge_descriptor h : curve_halfedges)
      {
        put(visited_edges, edge(h, mesh_), true);

        curve_length += CGAL::sqrt(CGAL::squared_distance(get(vpmap_, source(h, mesh_)),
                                                          get(vpmap_, target(h, mesh_))));

      }

      std::cout << "curve of length: " << curve_length << " (bound: " << feature_length
                << ") made of " << curve_halfedges.size() << " edges" << std::endl;

      if(curve_length <= feature_length)
        untag_curve(curve_halfedges);
    }

    // Now, loops are left
    for(halfedge_descriptor start_h : halfedges(mesh_))
    {
      if(!get(candidate_edges_, edge(start_h, mesh_)) || get(visited_edges, edge(start_h, mesh_)))
        continue;

      std::cout << "Start a new cycle at " << edge(start_h, mesh_) << std::endl;

      // simplified version of 'walk_curve' since we know these are cycles
      std::list<halfedge_descriptor> curve_halfedges;

      halfedge_descriptor curr_h = start_h;
      do
      {
        curve_halfedges.push_back(curr_h);
        curr_h = next_candidate_edge(curr_h);
      }
      while (curr_h != start_h);

      double curve_length = 0.;
      for(halfedge_descriptor h : curve_halfedges)
      {
        put(visited_edges, edge(curr_h, mesh_), true);

        curve_length += CGAL::sqrt(CGAL::squared_distance(get(vpmap_, source(h, mesh_)),
                                                          get(vpmap_, target(h, mesh_))));
      }

      std::cout << "cycle of length: " << curve_length << " (bound: " << feature_length
                << ") made of " << curve_halfedges.size() << " edges" << std::endl;

      if(curve_length <= feature_length)
        untag_curve(curve_halfedges);
    }
  }

  void initialize_angle_bounds(const FT theta_F, const FT theta_f, const FT theta_D, const FT theta_T,
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
  void tag_sharp_edges(const double feature_length,
                       const bool use_upper_bound_on_DA,
                       // all angles in degrees
                       const FT theta_F,
                       const FT theta_f,
                       const FT theta_D = 60,
                       const FT theta_T = 120,
                       const FT theta_t = 20,
                       const FT theta_e = 25,
                       const FT theta_k = 50,
                       int k = 5)
  {
    std::cout << "Tagging sharp edges..." << std::endl;

    use_upper_bound_on_DA_ = use_upper_bound_on_DA;

    initialize_angle_bounds(theta_F, theta_f, theta_D, theta_T, theta_t, theta_e, theta_k, k);

    initialize_maps();

    tag_candidate_edges(); // also builds a list of incident candidate halfedges at each vertex

    trim_candidate_list(feature_length);

    // In a post-phase, get rid of sharp polylines with a total length that is small
    remove_small_features(feature_length);

    // Mark remaining candidates halfedges as C1 discontinuities
#ifndef SHARP_EDGES_IN_INDEPENDENT_FILES
    std::ofstream out("results/all_sharp_edges.polylines.txt");
#endif
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

  double feature_size_;

  bool use_upper_bound_on_DA_;
};

template <typename TriangleMesh, typename VPM, typename GeomTraits>
double default_feature_length(const TriangleMesh& tmesh,
                              const VPM vpmap,
                              const GeomTraits& gt)
{
  const Bbox_3 bb = bbox(tmesh, CGAL::parameters::geom_traits(gt).vertex_point_map(vpmap));

  const double bbox_diagonal = CGAL::sqrt(CGAL::square(bb.xmax() - bb.xmin()) +
                                          CGAL::square(bb.ymax() - bb.ymin()) +
                                          CGAL::square(bb.zmax() - bb.zmin()));

  const double def_val = 0.01 * bbox_diagonal; // default filter is 1% of the bbox's diagonal

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
  std::cout << "default feature length: " << def_val << std::endl;
#endif

  return def_val;
}

} // end namespace internal

namespace experimental {

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
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type                VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, pmesh));

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type                    GeomTraits;
  GeomTraits gt = choose_parameter(get_parameter(np, internal_np::geom_traits), GeomTraits());

  typedef typename internal::Detector<PolygonMesh, VPM, GeomTraits, EdgeIsFeatureMap>   Detector;
  Detector detector(pmesh, vpm, gt, edge_is_feature_map);

  const double feature_length = internal::default_feature_length(pmesh, vpm, gt);

  const double weak_DA_in_deg = choose_parameter(get_parameter(np, internal_np::weak_dihedral_angle), 10.);
  const bool use_upper_bound_on_DA = choose_parameter(get_parameter(np, internal_np::use_upper_DA_bound), false);

#ifdef CGAL_PMP_DETECT_FEATURES_PP_DEBUG
  std::cout << "Use upper bound on DA? " << use_upper_bound_on_DA << std::endl;
#endif

  detector.tag_sharp_edges(feature_length, use_upper_bound_on_DA, strong_DA_in_deg, weak_DA_in_deg);
}

template <typename PolygonMesh,
          typename FT,
          typename EdgeIsFeatureMap,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
void detect_sharp_edges_pp(PolygonMesh& pmesh,
                           const FT strong_DA_in_deg,
                           EdgeIsFeatureMap edge_is_feature_map,
                           const CGAL_NP_CLASS& np)
{
  return detect_sharp_edges_pp(faces(pmesh), pmesh, strong_DA_in_deg, edge_is_feature_map, np);
}

template <typename PolygonMesh,
          typename FT,
          typename EdgeIsFeatureMap>
void detect_sharp_edges_pp(PolygonMesh& pmesh,
                           const FT strong_DA_in_deg,
                           EdgeIsFeatureMap edge_is_feature_map)
{
  return detect_sharp_edges_pp(faces(pmesh), pmesh, strong_DA_in_deg, edge_is_feature_map,
                               CGAL::parameters::all_default());
}

} // end namespace experimental
} // end namespace Polygon_mesh_processing
} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_DETECT_FEATURES_IN_POLYGON_MESH_PP_H
