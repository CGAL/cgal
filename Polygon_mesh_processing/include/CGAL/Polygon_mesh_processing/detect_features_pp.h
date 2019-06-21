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

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

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

  typedef CGAL::dynamic_vertex_property_t<double>                                 Vertex_double_tag;
  typedef CGAL::dynamic_vertex_property_t<bool>                                   Vertex_bool_tag;
  typedef CGAL::dynamic_halfedge_property_t<double>                               Halfedge_double_tag;
  typedef CGAL::dynamic_edge_property_t<double>                                   Edge_double_tag;
  typedef CGAL::dynamic_edge_property_t<bool>                                     Edge_bool_tag;

  typedef typename boost::property_map<PolygonMesh, Vertex_double_tag>::type      Vertex_double_pmap;
  typedef typename boost::property_map<PolygonMesh, Vertex_bool_tag>::type        Vertex_bool_pmap;
  typedef typename boost::property_map<PolygonMesh, Halfedge_double_tag>::type    Halfedge_double_pmap;
  typedef typename boost::property_map<PolygonMesh, Edge_double_tag>::type        Edge_double_pmap;
  typedef typename boost::property_map<PolygonMesh, Edge_bool_tag>::type          Edge_bool_pmap;

  typedef CGAL::Eigen_diagonalize_traits<double, 3>                               Diagonalize_traits;
  typedef Eigen::Matrix<double, 3, 3>                                             EigenMatrix;
  typedef Eigen::Matrix<double, 3, 1>                                             EigenVector;

public:
  Detector(PolygonMesh& mesh,
           const VPM vpmap,
           const GeomTraits& traits,
           const EIFMap ifemap /*is_feature_edge_map*/)
    : mesh_(mesh), vpmap_(vpmap), gt_(traits), ifemap_(ifemap)
  {
    CGAL_assertion(theta_T_ >= 2 * theta_t_);
    CGAL_assertion(theta_f_ <= theta_e_ && theta_e_ <= theta_F_);

    angle_defects_ = get(Vertex_double_tag(), mesh_);
    turning_angles_ = get(Vertex_double_tag(), mesh_);
    osta_angles_ = get(Halfedge_double_tag(), mesh_);
    dihedral_angles_ = get(Edge_double_tag(), mesh_);

    candidate_valence_ = get(Vertex_bool_tag(), mesh_);
    candidate_edges_ = get(Edge_bool_tag(), mesh_);
  }

private:
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
      return 2 * (CGAL_PI - total_angle);
    else
      return 2 * CGAL_PI - total_angle;
  }

  double dihedral_angle(const halfedge_descriptor h) const
  {
    if(is_border_edge(h, mesh_))
      return CGAL_PI;

    return gt_.compute_approximate_dihedral_angle_3_object()(get(vpmap_, source(h, mesh_)),
                                                             get(vpmap_, target(h, mesh_)),
                                                             get(vpmap_, target(next(h, mesh_), mesh_)),
                                                             get(vpmap_, target(next(opposite(h, mesh_), mesh_), mesh_)));
  }

  double dihedral_angle(const edge_descriptor e) const { return dihedral_angle(halfedge(e, mesh_)); }

  double turning_angle(const halfedge_descriptor first_h, const halfedge_descriptor second_h) const
  {
    return angle(get(vpmap_, target(first_h, mesh_)),
                 get(vpmap_, source(first_h, mesh_)),
                 get(vpmap_, target(second_h, mesh_)));
  }

  // obviously @cache that
  Vector ridge_direction(const vertex_descriptor v) const
  {
    // construct the matrix
    EigenMatrix M;

    for(halfedge_descriptor h : CGAL::halfedges_around_target(v, mesh_))
    {
      if(is_border(h, mesh_))
        continue;

      const face_descriptor f = face(h, mesh_);

      const Vector fn = Polygon_mesh_processing::compute_face_normal(f, mesh_, CGAL::parameters::vertex_point_map(vpmap_)
                                                                                                .geom_traits(gt_));

      const double face_area = Polygon_mesh_processing::face_area(f, mesh_, CGAL::parameters::vertex_point_map(vpmap_)
                                                                                             .geom_traits(gt_));

      EigenVector mj;
      mj << fn[0] << fn[1] << fn[2];

      M += face_area * mj * mj.transpose();
    }

    // compute the eigenvalues/vectors
    Eigen::SelfAdjointEigenSolver<EigenMatrix> eigen_solver;

    eigen_solver.computeDirect(M);

    // sort eigenvalues (& vectors)

    // check degeneracy

    // get ridge direction

  }

  double one_side_turning_angle(const halfedge_descriptor h) const
  {
    const vertex_descriptor v = source(h, mesh_);
    const Vector dir = ridge_direction(v);
    const Point r = gt_.construct_translated_point_3_object()(get(vpmap_, v), dir);

    return turning_angle(get(vpmap_, target(h, mesh_)), get(vpmap_, v), r);
  }

  // -----------------------------------------------------------------------------------------------
  // ------------------------------------- OSTA/DA STRENGTH ----------------------------------------
  // -----------------------------------------------------------------------------------------------
  bool is_locally_strong_in_DA(const halfedge_descriptor h) const
  {
    if(is_border_edge(h, mesh_))
      return true;

    const double dih_angle = get(dihedral_angles_, edge(h, mesh_));
    if(dih_angle < theta_f_)
      return false;

    double max_other_DA_osta_smaller_than_half_pi = - std::numeric_limits<double>::max();
    double max_other_DA_osta_greater_than_half_pi = - std::numeric_limits<double>::max();
    for(halfedge_descriptor other_h : CGAL::halfedges_around_source(h, mesh_))
    {
      if(other_h == h)
        continue;

      const edge_descriptor other_e = edge(other_h, mesh_);
      const double other_dih_angle = get(dihedral_angles_, other_e);
      if(other_dih_angle < theta_f_)
        continue;

      const double twice_other_osta = 2 * get(osta_angles_, other_h);

      // @fixme really ignore osta == pi/2 ?
      if(twice_other_osta < CGAL_PI) // <=> osta < pi/2
        max_other_DA_osta_smaller_than_half_pi = (CGAL::max)(max_other_DA_osta_smaller_than_half_pi,
                                                             other_dih_angle);
      else if(twice_other_osta > CGAL_PI) // <=> osta < pi/2
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

    if(osta_angle < theta_t_)
    {
      double min_osta_value = std::numeric_limits<double>::max();
      for(halfedge_descriptor other_h : CGAL::halfedges_around_source(h, mesh_)) // @fixme? border
      {
        if(other_h == h)
          continue;

        const double other_DA = get(dihedral_angles_, edge(other_h, mesh_));
        if(other_DA < theta_f_)
          continue;

        min_osta_value = (CGAL::min)(min_osta_value, get(osta_angles_, h));
      }

      if(osta_angle <= min_osta_value)
        return true;
    }
    else if(osta_angle > CGAL_PI - theta_t_) // @strict
    {
      double max_osta_value = -(std::numeric_limits<double>::max());
      for(halfedge_descriptor other_h : CGAL::halfedges_around_source(h, mesh_)) // @fixme? border
      {
        if(other_h == h)
          continue;

        const double other_DA = get(dihedral_angles_, edge(other_h, mesh_));
        if(other_DA < theta_f_)
          continue;

        max_osta_value = (CGAL::max)(max_osta_value, get(osta_angles_, h));
      }

      if(osta_angle >= max_osta_value)
        return true;
    }

    return false;
  }

  bool is_strong_in_DA(const edge_descriptor e) const { return get(dihedral_angles_, e) > theta_e_; } // e-strong
  bool is_unconditionally_strong_in_DA(const vertex_descriptor v) const { return get(angle_defects_, v) > theta_D_; }
  bool is_unconditionally_strong_in_DA(const edge_descriptor e) const { return get(dihedral_angles_, e) > theta_F_; }
  bool is_unconditionally_strong_in_TA(const vertex_descriptor v) const { return get(turning_angles_, v) > theta_T_; }

  bool is_ambiguous_vertex(const vertex_descriptor v) const { } // @todo
  bool is_sharp_corner(const vertex_descriptor v) const { return is_unconditionally_strong_in_DA(v); }
  bool is_sharp_edge(const edge_descriptor e) const { return is_unconditionally_strong_in_DA(e); }

  // -----------------------------------------------------------------------------------------------
  // --------------------------------------- ATTACHEDNESS ------------------------------------------
  // -----------------------------------------------------------------------------------------------
  bool is_attached(const halfedge_descriptor h) const
  {
    const double dih_angle = get(dihedral_angles_, edge(h, mesh_));

    // @fixme weird formulation in the paper (p. 13 bottom) OR OR OR? is the (1) misplaced?
    if(dih_angle < theta_f_)
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
    const double dih_angle = get(dihedral_angles_, edge(h, mesh_));

    if(dih_angle < theta_f_)
      return false;

    const vertex_descriptor v = source(h, mesh_);
    if(is_sharp_edge(edge(h, mesh_)) || is_sharp_corner(v) || is_ambiguous_vertex(v))
      return true;

    if(is_locally_strong_in_OSTA(h) && is_locally_strong_in_DA(h)) // note the '&&'
      return true;
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
    return (get(candidate_edges_, edge(h, mesh_)) &&
            get(candidate_valence_, source(h, mesh_)) == 1);
  }

  bool is_dangling(const halfedge_descriptor h) const
  {
    if(!is_singleton(h))
      return false;

    return (!is_unconditionally_strong_in_DA(edge(h, mesh_)) ||
            !is_sharp_corner(source(h, mesh_)));
  }

  bool is_semi_joint(const halfedge_descriptor h) const
  {
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
          second_h = other_h;
        else
          return false; // strictly more than 2 incident candidate halfedges
      }
    }

    if(second_h == boost::graph_traits<PolygonMesh>::null_halfedge()) // only 1 incident candidate edge
      return false;

    // Now we have exactly 2 candidate edges
    edge_descriptor first_e = edge(first_h, mesh_);
    edge_descriptor second_e = edge(second_h, mesh_);

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

      if(2 * get(dihedral_angles_, edge(h, mesh_)) > CGAL_PI) // @todo >= (check across the whole file actually)
        return true;
    }

    return false;
  }

  bool is_disjoint(const halfedge_descriptor h) const
  {
    const vertex_descriptor v = source(h, mesh_);
    if(get(candidate_valence_, v) < 3)
      return false;

    // @cache stuff (esp. 'is_locally_strong_in_...')

    const bool cond_1 = (!is_locally_strong_in_DA(h) &&
                         !is_locally_strong_in_OSTA(h) &&
                         !is_unconditionally_strong_in_DA(edge(h, mesh_)) &&
                         !is_sharp_corner(v) &&
                         !is_ambiguous_vertex(v));

    const bool cond_2 = (!is_locally_strong_in_DA(h) ||
                         !is_locally_strong_in_OSTA(h) &&
                         !is_strong_in_DA(edge(h, mesh_)) &&
                         !is_sharp_corner(v) &&
                         !is_ambiguous_vertex(v));

    const bool cond_3 = (!is_strong_in_DA(edge(h, mesh_)) &&
                         has_incident_unconditionally_strong_in_DA_incident_edges(v));

    const bool cond_4 = (!is_unconditionally_strong_in_DA(edge(h, mesh_)) &&
                         has_incident_acute_edge(v));
  }

  bool is_multi_joint(const halfedge_descriptor h) const
  {
    return (get(candidate_valence_, source(h, mesh_)) >= 3 && !is_disjoint(h));
  }

  bool is_end_halfedge(const halfedge_descriptor h) const
  {
    // @cache all that stuff
    return (is_singleton(h) || is_semi_joint(h) || is_disjoint(h) || is_multi_joint(h));
  }

  bool is_obscure_end_halfedge(const halfedge_descriptor h) const
  {
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
    if(is_salient_curve(start, last))
      return false;

    const halfedge_descriptor last_opp = opposite(last, mesh_);
    const bool is_start_end_obscure = is_obscure_end_halfedge(start);
    const bool is_last_end_obscure = is_obscure_end_halfedge(last_opp);

    const bool cond_1 = ((is_start_end_obscure && is_last_end_obscure) ||
                         is_dangling(start) || is_dangling(last_opp)) && // @fixme xor?
                        number_of_edges_with_DAs_over_threshold < k_; // @fixme, with equality?
    if(cond_1)
      return true;

    const bool cond_2 = (is_start_end_obscure || is_last_end_obscure) && // @fixme xor?
                        has_incident_unconditionally_strong_in_DA_incident_edges(source(start, mesh_)) &&
                        has_incident_unconditionally_strong_in_DA_incident_edges(target(last, mesh_)) &&
                        number_of_edges_with_DAs_over_threshold == 0;

    return cond_2;
  }

  bool is_semi_salient_curve(const halfedge_descriptor start, const halfedge_descriptor last) const
  {
    return !is_salient_curve(start, last) && !is_obscure_curve(start, last);
  }

  void tag_candidate_edges()
  {
    for(vertex_descriptor v : vertices(mesh_))
      put(candidate_valence_, v, 0);

    for(edge_descriptor e : edges(mesh_))
    {
      const bool is_qs_edge = is_quasi_strong(e);
      put(candidate_edges_, e, is_qs_edge);

      if(is_qs_edge)
      {
        put(candidate_valence_, source(e, mesh_), get(candidate_valence_, source(e, mesh_)) + 1);
        put(candidate_valence_, target(e, mesh_), get(candidate_valence_, target(e, mesh_)) + 1);
      }
    }
  }

  // @fixme completely ignoring the fact that sharp corners are not marked as candidates currently
  // (don't really care since we want to tag sharp edges, but might want to be able to tag
  // sharp vertices too).

  halfedge_descriptor next_candidate_edge(const halfedge_descriptor in_h) const
  {
    CGAL_precondition(get(candidate_valence_, target(in_h, mesh_)) == 2);

    const edge_descriptor in_e = edge(in_h, mesh_);
    for(halfedge_descriptor h : CGAL::halfedges_around_source(target(in_h, mesh_), mesh_))
    {
      if(edge(h, mesh_) == in_e)
        continue;

      if(get(candidate_edges_, edge(h, mesh_)))
        return h;
    }

    CGAL_assertion(false);
    return boost::graph_traits<PolygonMesh>::null_halfedge();
  }

  void walk_curve(const halfedge_descriptor start,
                  std::list<halfedge_descriptor>& curve_halfedges,
                  int& number_of_edges_with_DAs_over_threshold) const
  {
    halfedge_descriptor curr_h = start;
    for(;;)
    {
      CGAL_assertion(get(candidate_edges_, edge(curr_h, mesh_)));
      curve_halfedges.push_back(curr_h);

      if(get(dihedral_angles_, edge(curr_h, mesh_)) > theta_k_)
        ++number_of_edges_with_DAs_over_threshold;

      // is this correct?
      if(is_end_halfedge(opposite(curr_h, mesh_)))
        break;

      CGAL_assertion_code(halfedge_descriptor old_h = curr_h;)
      curr_h = next_candidate_edge(curr_h);
      CGAL_assertion(old_h != curr_h);
      CGAL_assertion(target(old_h, mesh_) == source(curr_h, mesh_));

      if(is_end_halfedge(curr_h)) // meaning the next one is _not_ added to the walk
        break;
    }
  }

  void initialize_maps()
  {
    for(vertex_descriptor v : vertices(mesh_))
    {
      put(angle_defects_, v, angle_defect(v));
      const Vector dir = ridge_direction();
    }

    for(halfedge_descriptor h : halfedges(mesh_))
      put(dihedral_angles_, h, dihedral_angle(h));

    for(halfedge_descriptor h : halfedges(mesh_))
    {
      put(osta_angles_, h, one_side_turning_angle(h));

      // @cache local extremas at vertices (don't forget about the DA filter)
    }
  }

public:
  bool tag_sharp_edges()
  {
    // build list of incident candidate halfedges at each vertex
    tag_candidate_edges();

    bool removed_some_halfedges = false;
    do
    {
      // collect obscure end-edges from quasi-strong edges
      std::unordered_set<halfedge_descriptor> obscure_end_halfedges;
      for(halfedge_descriptor h : halfedges(mesh_))
      {
        // @cache order that so it's cheap
        if(get(candidate_edges_, edge(h, mesh_)) && is_obscure_end_halfedge(h))
          obscure_end_halfedges.insert(h);
      }

      // while there are still obscure end-edges left
      while(!obscure_end_halfedges.empty())
      {
        // traverse curve
        halfedge_descriptor start = *(obscure_end_halfedges.begin());

        std::list<halfedge_descriptor> curve_halfedges;
        int number_of_edges_with_DAs_over_threshold = 0;

        walk_curve(start, curve_halfedges, number_of_edges_with_DAs_over_threshold);

        const halfedge_descriptor last_h = curve_halfedges.back();
        // check if the curve is obscure
        if(is_obscure_curve(start, last_h, number_of_edges_with_DAs_over_threshold))
        {
          halfedge_descriptor h = start;
          for(;;)
          {
            put(candidate_edges_, edge(h, mesh_), false);
            CGAL_assertion(get(candidate_valence_, source(h, mesh_)) > 0);
            put(candidate_valence_, source(h, mesh_), get(candidate_valence_, source(h, mesh_)) - 1);

            if(h == last_h)
              break;
          }

          CGAL_assertion(get(candidate_valence_, target(last_h, mesh_)) > 0);
          put(candidate_valence_, target(last_h, mesh_), get(candidate_valence_, target(last_h, mesh_)) - 1);

          // Reasoning is that, tautologically, we cannot have met any obscure end halfedge
          // before, otherwise we would have stopped
          obscure_end_halfedges.erase(opposite(last_h, mesh_));

          removed_some_halfedges = true;
        }
      }
    }
    while(removed_some_halfedges); // while removed halfedges

    // mark remaining candidates halfedges as C1 discontinuities
  }

private:
  PolygonMesh& mesh_;
  const VPM vpmap_;
  GeomTraits gt_;
  const EIFMap ifemap_;

  // @todo all of that into dynamic pmaps and 'at' --> get()
  Vertex_double_pmap angle_defects_;
  Vertex_double_pmap turning_angles_;
  Halfedge_double_pmap osta_angles_;
  Edge_double_pmap dihedral_angles_;

  Vertex_bool_pmap candidate_valence_;
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

template <typename PolygonMesh,
          typename FaceRange,
          typename EdgeIsFeatureMap,
          typename NamedParameters>
void detect_sharp_edges_pp(const FaceRange& /*faces*/,
                           PolygonMesh& pmesh,
                           EdgeIsFeatureMap edge_is_feature_map,
                           const NamedParameters& np)
{
  using boost::choose_param;
  using boost::get_param;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type          VPM;
  VPM vpm = choose_param(get_param(np, internal_np::vertex_point),
                         get_const_property_map(CGAL::vertex_point, pmesh));

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type              GeomTraits;
  GeomTraits traits = choose_param(get_param(np, internal_np::geom_traits), GeomTraits());

  internal::Detector<PolygonMesh, VPM, GeomTraits, EdgeIsFeatureMap> detector(pmesh, vpm, traits, edge_is_feature_map);

  detector.tag_sharp_edges();
}

template <typename PolygonMesh,
          typename EdgeIsFeatureMap,
          typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
void detect_sharp_edges_pp(PolygonMesh& pmesh,
                           EdgeIsFeatureMap edge_is_feature_map,
                           const CGAL_PMP_NP_CLASS& np)
{
  return detect_sharp_edges_pp(faces(pmesh), pmesh, edge_is_feature_map, np);
}

} // end namespace PMP
} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYGON_MESH_PP_H
