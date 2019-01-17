// Copyright (c) 2019  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Kaimo Hu

#ifndef REMESH_H_
#define REMESH_H_

// STL
#include <map>
// local
#include "types.h"
#include "pvertex.h"
#include "phalfedge.h"
#include "pfacet.h"

template <class Kernel, class Polyhedron>
class Remesh {
 public:
  // types
  typedef typename Polyhedron::Vector Vector;
  typedef typename Polyhedron::Normal Normal;
  // iterators
  typedef typename Polyhedron::Edge_iterator Edge_iterator;
  typedef typename Polyhedron::Vertex_iterator Vertex_iterator;
  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  typedef typename Polyhedron::Halfedge_iterator Halfedge_iterator;
  // handles
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Polyhedron::Vertex_handle Vertex_handle;
  typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;
  // circulators
  typedef typename Polyhedron::Halfedge_around_vertex_circulator
                               Halfedge_around_vertex_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_const_circulator
                               Halfedge_around_vertex_const_circulator;
  // point list
  typedef typename Polyhedron::Point_Comp Point_Comp;
  typedef typename Polyhedron::Point Point;
  typedef typename Polyhedron::Point_list Point_list;
  typedef typename Polyhedron::Point_iter Point_iter;
  typedef typename Polyhedron::Point_const_iter Point_const_iter;
  typedef typename Polyhedron::Halfedge_list Halfedge_list;
  typedef typename Polyhedron::Halfedge_iter Halfedge_iter;
  typedef typename Polyhedron::Halfedge_const_list Halfedge_const_list;
  typedef typename Polyhedron::Halfedge_const_iter Halfedge_const_iter;
  typedef typename Polyhedron::Facet_list Facet_list;
  typedef typename Polyhedron::Facet_iter Facet_iter;
  typedef typename Polyhedron::Facet_const_list Facet_const_list;
  typedef typename Polyhedron::Facet_const_iter Facet_const_iter;
  // Dynamic priority queues
  typedef CPHalfedge<FT, Halfedge_handle> PHalfedge;
  typedef CDPQueue_long<PHalfedge> DPQueue_halfedge_long;
  typedef CDPQueue_short<PHalfedge> DPQueue_halfedge_short;
  typedef CPVertex<FT, Vertex_handle> PVertex;
  typedef CDPQueue_long<PVertex> DPQueue_vertex_long;
  typedef CDPQueue_short<PVertex> DPQueue_vertex_short;
  typedef CPFacet<FT, Facet_handle> PFacet;
  typedef CDPQueue_long<PFacet> DPQueue_facet_long;
  typedef CDPQueue_short<PFacet> DPQueue_facet_short;
  // visit list and iterator
  typedef typename std::list<std::pair<Point, FT>> Visit_list;
  typedef typename std::list<std::pair<Point, FT>>::iterator Visit_iter;

  // life cycle
  Remesh() {
    m_max_error_threshold = 0.2;          // general paramters
    m_min_angle_threshold = 30.0;
    m_max_mesh_complexity = 100000000;
    m_smooth_angle_delta = 0.1;
    m_apply_edge_flip = true;
    m_edge_flip_strategy = EdgeFlipStrategy::k_improve_angle;
    m_flip_after_split_and_collapse = true;
    m_relocate_after_local_operations = true;
    m_relocate_strategy = RelocateStrategy::k_cvt_barycenter;
    m_keep_vertex_in_one_ring = false;
    m_use_local_aabb_tree = true;
    m_collapsed_list_size = 10;
    m_decrease_max_errors = true;
    m_track_information = true;
    m_apply_initial_mesh_simplification = true;
    m_apply_final_vertex_relocation = true;

    m_samples_per_facet_in = 10;          // sample parameters
    m_samples_per_facet_out = 10;
    m_max_samples_per_area = 10000;
    m_min_samples_per_triangle = 1;
    m_bvd_iteration_count = 1;
    m_sample_number_strategy = SampleNumberStrategy::k_fixed;
    m_sample_strategy = SampleStrategy::k_adaptive;
    m_use_stratified_sampling = false;

    m_sum_theta = 1.0;                    // feature parameters
    m_sum_delta = 0.5;
    m_dihedral_theta = 1.0;
    m_dihedral_delta = 0.5;
    m_feature_difference_delta = 0.15;
    m_feature_control_delta = 0.5;
    m_inherit_element_types = false;
    m_use_feature_intensity_weights = false;

    m_vertex_optimize_count = 2;          // vertex optimization parameters
    m_vertex_optimize_ratio = 0.9;
    m_stencil_ring_size = 1;
    m_optimize_strategy = OptimizeStrategy::k_approximation;
    m_facet_optimize_type = OptimizeType::k_both;
    m_edge_optimize_type = OptimizeType::k_both;
    m_vertex_optimize_type = OptimizeType::k_both;
    m_optimize_after_local_operations = true;

    initialize_private_data();
  }

  Remesh(FT max_error_threshold,  // general parameters
    FT min_angle_threshold,
    int max_mesh_complexity,
    FT smooth_angle_delta,
    bool apply_edge_flip,
    EdgeFlipStrategy edge_flip_strategy,
    bool flip_after_split_and_collapse,
    bool relocate_after_local_operations,
    RelocateStrategy relocate_strategy,
    bool keep_vertex_in_one_ring,
    bool use_local_aabb_tree,
    int collapse_list_size,
    bool decrease_max_errors,
    bool track_information,
    bool apply_initial_mesh_simplification,
    bool apply_final_vertex_relocation,
    int samples_per_facet_in,         // sample parameters
    int samples_per_facet_out,
    int max_samples_per_area,
    int min_samples_per_triangle,
    int bvd_iteration_count,
    SampleNumberStrategy sample_number_strategy,
    SampleStrategy sample_strategy,
    bool use_stratified_sampling,
    FT sum_theta,                 // feature intensity parameters
    FT sum_delta,
    FT dihedral_theta,
    FT dihedral_delta,
    FT feature_difference_delta,
    FT feature_control_delta,
    bool inherit_element_types,
    bool use_feature_intensity_weights,
    int vertex_optimize_count,      // vertex optimization parameters
    FT vertex_optimize_ratio,
    int stencil_ring_size,
    OptimizeStrategy optimize_strategy,
    OptimizeType facet_optimize_type,
    OptimizeType edge_optimize_type,
    OptimizeType vertex_optimize_type,
    bool optimize_after_local_operations)
    : m_max_error_threshold(max_error_threshold),
    m_min_angle_threshold(min_angle_threshold),
    m_max_mesh_complexity(max_mesh_complexity),
    m_smooth_angle_delta(smooth_angle_delta),
    m_apply_edge_flip(apply_edge_flip),
    m_edge_flip_strategy(edge_flip_strategy),
    m_flip_after_split_and_collapse(flip_after_split_and_collapse),
    m_relocate_after_local_operations(relocate_after_local_operations),
    m_relocate_strategy(relocate_strategy),
    m_keep_vertex_in_one_ring(keep_vertex_in_one_ring),
    m_use_local_aabb_tree(use_local_aabb_tree),
//    m_collapse_list_size(collapse_list_size),
    m_decrease_max_errors(decrease_max_errors),
    m_track_information(track_information),
//    m_apply_initial_mesh_simplificaiton(apply_initial_mesh_simplification),
    m_apply_final_vertex_relocation(apply_final_vertex_relocation),
    m_samples_per_facet_in(samples_per_facet_in),
    m_samples_per_facet_out(samples_per_facet_out),
    m_max_samples_per_area(max_samples_per_area),
    m_min_samples_per_triangle(min_samples_per_triangle),
    m_bvd_iteration_count(bvd_iteration_count),
    m_sample_number_strategy(sample_number_strategy),
    m_sample_strategy(sample_strategy),
    m_use_stratified_sampling(use_stratified_sampling),
    m_sum_theta(sum_theta),
    m_sum_delta(sum_delta),
    m_dihedral_theta(dihedral_theta),
    m_dihedral_delta(dihedral_delta),
    m_feature_difference_delta(feature_difference_delta),
    m_feature_control_delta(feature_control_delta),
    m_inherit_element_types(inherit_element_types),
    m_use_feature_intensity_weights(use_feature_intensity_weights),
    m_vertex_optimize_count(vertex_optimize_count),
    m_vertex_optimize_ratio(vertex_optimize_ratio),
    m_stencil_ring_size(stencil_ring_size),
    m_optimize_strategy(optimize_strategy),
    m_facet_optimize_type(facet_optimize_type),
    m_edge_optimize_type(edge_optimize_type),
    m_vertex_optimize_type(vertex_optimize_type),
    m_optimize_after_local_operations(optimize_after_local_operations) {

    initialize_private_data();
  }

  ~Remesh() {}

 public:
  // ) parameter access functions
  FT get_max_error_threshold() const { return m_max_error_threshold; }
  void set_max_error_threshold(FT value) { m_max_error_threshold = value; }
  FT get_min_angle_threshold() const { return m_min_angle_threshold; }
  void set_min_angle_threshold(FT value) { m_min_angle_threshold = value; }
  int get_max_mesh_complexity() const { return m_max_mesh_complexity; }
  void set_max_mesh_complexity(int value) { m_max_mesh_complexity = value; }
  FT get_smooth_angle_delta() const { return m_smooth_angle_delta; }
  void set_smooth_angle_delta(FT value) { m_smooth_angle_delta = value; }
  bool get_apply_edge_flip() const { return m_apply_edge_flip; }
  void set_apply_edge_flip(bool value) { m_apply_edge_flip = value; }
  EdgeFlipStrategy get_edge_flip_strategy() const
      { return m_edge_flip_strategy; }
  void set_edge_flip_strategy(EdgeFlipStrategy value)
      { m_edge_flip_strategy = value; }
  bool get_flip_after_split_and_collapse() const
      { return m_flip_after_split_and_collapse; }
  void set_flip_after_split_and_collapse(bool value)
      { m_flip_after_split_and_collapse = value; }
  bool get_relocate_after_local_operations() const
      { return m_relocate_after_local_operations; }
  void set_relocate_after_local_operations(bool value)
      { m_relocate_after_local_operations = value; }
  RelocateStrategy get_relocate_strategy() const
      { return m_relocate_strategy; }
  void set_relocate_strategy(RelocateStrategy value)
      { m_relocate_strategy = value; }
  bool get_keep_vertex_in_one_ring() const
      { return m_keep_vertex_in_one_ring; }
  void set_keep_vertex_in_one_ring(bool value)
      { m_keep_vertex_in_one_ring = value; }
  bool get_use_local_aabb_tree() const { return m_use_local_aabb_tree; }
  void set_use_local_aabb_tree(bool value) { m_use_local_aabb_tree = value; }
  int get_collapsed_list_size() const { return m_collapsed_list_size; }
  void set_collapsed_list_size(int value) { m_collapsed_list_size = value; }
  bool get_decrease_max_errors() const { return m_decrease_max_errors; }
  void set_decrease_max_errors(bool value) { m_decrease_max_errors = value; }
  bool get_track_information() const { return m_track_information; }
  void set_track_information(bool value) { m_track_information = value; }
  bool get_apply_initial_mesh_simplification() const
      { return m_apply_initial_mesh_simplification; }
  void set_apply_initial_mesh_simplification(bool value)
      { m_apply_initial_mesh_simplification = value; }
  bool get_apply_final_vertex_relocation() const
      { return m_apply_final_vertex_relocation; }
  void set_apply_final_vertex_relocation(bool value)
      { m_apply_final_vertex_relocation = value; }

  int get_samples_per_facet_in() const { return m_samples_per_facet_in; }
  void set_samples_per_facet_in(int value) { m_samples_per_facet_in = value; }
  int get_samples_per_facet_out() const { return m_samples_per_facet_out; }
  void set_samples_per_facet_out(int value)
      { m_samples_per_facet_out = value; }
  int get_max_samples_per_area() const { return m_max_samples_per_area; }
  void set_max_samples_per_area(int value) { m_max_samples_per_area = value; }
  int get_min_samples_per_triangle() const
      { return m_min_samples_per_triangle; }
  void set_min_samples_per_triangle(int value)
      { m_min_samples_per_triangle = value; }
  int get_bvd_iteration_count() const { return m_bvd_iteration_count; }
  void set_bvd_iteration_count(int value) { m_bvd_iteration_count = value; }
  SampleNumberStrategy get_sample_number_strategy() const
      { return m_sample_number_strategy; }
  void set_sample_number_strategy(SampleNumberStrategy value)
      { m_sample_number_strategy = value; }
  SampleStrategy get_sample_strategy() const { return m_sample_strategy; }
  void set_sample_strategy(SampleStrategy value) { m_sample_strategy = value; }
  bool get_use_stratified_sampling() const
      { return m_use_stratified_sampling; }
  void set_use_stratified_sampling(bool value)
      { m_use_stratified_sampling = value; }

  FT get_sum_theta() const { return m_sum_theta; }
  void set_sum_theta(FT value) { m_sum_theta = value; }
  FT get_sum_delta() const { return m_sum_delta; }
  void set_sum_delta(FT value) { m_sum_delta = value; }
  FT get_dihedral_theta() const { return m_dihedral_theta; }
  void set_dihedral_theta(FT value) { m_dihedral_theta = value; }
  FT get_dihedral_delta() const { return m_dihedral_delta; }
  void set_dihedral_delta(FT value) { m_dihedral_delta = value; }
  FT get_feature_difference_delta() const
      { return m_feature_difference_delta; }
  void set_feature_difference_delta(FT value)
      { m_feature_difference_delta = value; }
  FT get_feature_control_delta() const { return m_feature_control_delta; }
  void set_feature_control_delta(FT value) { m_feature_control_delta = value; }
  bool get_inherit_element_types() const { return m_inherit_element_types; }
  void set_inherit_element_types(bool value)
      { m_inherit_element_types = value; }
  bool get_use_feature_intensity_weights() const
      { return m_use_feature_intensity_weights; }
  void set_use_feature_intensity_weights(bool value)
      { m_use_feature_intensity_weights = value; }

  int get_vertex_optimize_count() const { return m_vertex_optimize_count; }
  void set_vertex_optimize_count(int value)
      { m_vertex_optimize_count = value; }
  FT get_vertex_optimize_ratio() const { return m_vertex_optimize_ratio; }
  void set_vertex_optimize_ratio(FT value) { m_vertex_optimize_ratio = value; }
  int get_stencil_ring_size() const { return m_stencil_ring_size; }
  void set_stencil_ring_size(int value) { m_stencil_ring_size = value; }
  OptimizeStrategy get_optimize_strategy() const
      { return m_optimize_strategy; }
  void set_optimize_strategy(OptimizeStrategy value)
      { m_optimize_strategy = value; }
  OptimizeType get_facet_optimize_type() const
      { return m_facet_optimize_type; }
  void set_facet_optimize_type(OptimizeType value)
      { m_facet_optimize_type = value; }
  OptimizeType get_edge_optimize_type() const { return m_edge_optimize_type; }
  void set_edge_optimize_type(OptimizeType value)
      { m_edge_optimize_type = value; }
  OptimizeType get_vertex_optimize_type() const
      { return m_vertex_optimize_type; }
  void set_vertex_optimize_type(OptimizeType value)
      { m_vertex_optimize_type = value; }
  bool get_optimize_after_local_operations() const
      { return m_optimize_after_local_operations; }
  void set_optimize_after_local_operations(bool value)
      { m_optimize_after_local_operations = value; }

  void initialize_private_data() {
    m_collapsed_list.clear();
    m_collapsed_map.clear();
  }

  // polyhedron manipulations
  int eliminate_degenerated_facets(Polyhedron *polyhedron) const {
    FT radian_threshold = MIN_VALUE;
    // step 1: fill all the degenerated facets in the dynamic priority queue
    DPQueue_facet_long degenerated_facets;
    fill_degenerated_facets_queue(radian_threshold,
                                  &degenerated_facets, polyhedron);
    // step 2: eliminate these degenerations one by one
    int nb_eliminations = static_cast<int>(degenerated_facets.size());
    FT longest_squared_length = 0.0;
    while (!degenerated_facets.empty()) {
      PFacet facet = degenerated_facets.top();
      degenerated_facets.pop();
      Facet_handle fh = facet.facet();
      Halfedge_handle shortest_hh = polyhedron->get_shortest_halfedge(fh);
      // case 1: the triangle is acute, so we only need to collpase
      if (polyhedron->squared_length(shortest_hh) < SQUARED_MIN_VALUE) {
        // remove facets before collapse
        remove_incident_facets(shortest_hh, &degenerated_facets, polyhedron);
        // collapse the edge
        Point new_point = polyhedron->midpoint(shortest_hh);
        Vertex_handle v_joined = polyhedron->collapse_short_edge(new_point,
                                                                 shortest_hh);
        if (v_joined == NULL) {
          return -1;
        }
        // add incident facets to the queue if it is still generated facets
        add_circulator_degenerate_facets(v_joined, radian_threshold,
                                         &degenerated_facets, polyhedron);
      }
      // case 2: the triangle is obtuse, so we split the longest edge first
      else {
        Halfedge_handle longest_hh = polyhedron->get_longest_halfedge(fh);
        longest_squared_length = polyhedron->squared_length(longest_hh);
        Point new_point = polyhedron->get_opposite_vertex(longest_hh)->point();
        // remove facets before split
        degenerated_facets.remove(PFacet(fh));
        if (!longest_hh->opposite()->is_border()) {
          degenerated_facets.remove(PFacet(longest_hh->opposite()->facet()));
        }
        // split the longest halfedge
        Halfedge_handle hnew = polyhedron->split_long_edge(new_point,
                                                           longest_hh);
        // add the new degenerated facets after split
        add_circulator_degenerate_facets(hnew->vertex(), radian_threshold,
                                         &degenerated_facets, polyhedron);
      }
    }
    return nb_eliminations;
  }

  int split_long_edges(Polyhedron *polyhedron) const {
    FT average_length = polyhedron->get_average_length();
    FT max_length = 4.0 / 3.0 * average_length;
    FT max_squared_length = max_length * max_length;
    // fill dynamic priority queue with long edges in front
    DPQueue_halfedge_long long_edges;
    fill_queue_with_long_edges(max_squared_length, &long_edges, polyhedron);
    int nb_split = 0;
    while (!long_edges.empty()) {
      // get a copy of the candidate edge
      PHalfedge edge = long_edges.top();
      long_edges.pop();
      Halfedge_handle hh = edge.halfedge();
      // remove the hh and its opposite from queue
      long_edges.remove(PHalfedge(hh));
      long_edges.remove(PHalfedge(hh->opposite()));
      // split the specified long edge
      Point new_point = polyhedron->midpoint(hh);
      Halfedge_handle hnew = polyhedron->split_long_edge(new_point, hh);
      // update the queue
      add_circular_long_edges(max_squared_length,
        hnew->vertex(), &long_edges, polyhedron);
      ++nb_split;
    }
    return nb_split;
  }

  // isotropic remeshing
  void isotropic_remeshing(const Facet_tree &input_facet_tree,
                           const Bbox &bbox, Polyhedron *m_pRemesh) {
    CGAL::Timer timer;
    timer.start();
    std::cout << "Isotropic remeshing..." << std::endl;
    if (m_apply_initial_mesh_simplification) {
      std::cout << std::endl;
      initial_mesh_simplification(input_facet_tree, bbox, m_pRemesh);
    }
    std::cout << std::endl;
    maximize_minimal_angle(input_facet_tree, bbox, m_pRemesh);
    if (m_apply_final_vertex_relocation) {
      std::cout << std::endl;
      final_vertex_relocation(input_facet_tree, bbox, m_pRemesh);
    }
    std::cout << std::endl;
    std::cout << "Done, (total time is " << timer.time() << " s)" << std::endl;
  }

  void initial_mesh_simplification(const Facet_tree &input_facet_tree,
                                   const Bbox &bbox, Polyhedron *m_pRemesh) {
    /* for isotropic case, fill priority queue with collapsible edges.
    The priority has two choices:
    1) the error before collapse: E_{before}
    2) edge_length * opposite_angle (current implementation) */
    FT max_error_threshold = get_max_error_threshold_value(bbox);
    CGAL::Timer timer;
    timer.start();
    std::cout << "Initial mesh simplification..." << std::endl;
    std::cout << "(max error threshold = " << max_error_threshold
              << ")" << std::endl;
    DPQueue_halfedge_long large_error_queue;
    DPQueue_halfedge_short collapse_candidate_queue;
    fill_collapse_candidate_edges(max_error_threshold, &large_error_queue,
                                  &collapse_candidate_queue, m_pRemesh);
    unsigned int index = 0, nb_operations = 0;
    Halfedge_handle max_error_halfedge;
    FT max_error = 0.0;
    while (m_pRemesh->size_of_vertices() < m_max_mesh_complexity &&
      (!large_error_queue.empty() || !collapse_candidate_queue.empty())) {
      while (m_pRemesh->size_of_vertices() < m_max_mesh_complexity &&
          !large_error_queue.empty()) {
        PHalfedge edge = large_error_queue.top();
        large_error_queue.pop();
        max_error_halfedge = edge.halfedge();
        max_error = CGAL::sqrt(edge.priority());
        if (m_track_information) {
          std::cout << ++index << ": error queue size = "
            << large_error_queue.size() << " ";
        }
        greedy_reduce_error(input_facet_tree, max_error_threshold, max_error,
            m_track_information, true, &large_error_queue,
            &collapse_candidate_queue, max_error_halfedge, m_pRemesh);
      }
      if (!collapse_candidate_queue.empty()) {
        if (m_track_information) {
          std::cout << ++index << ": collapse queue size = "
                    << collapse_candidate_queue.size() << " ";
        }
        // step 1: get the top halfedge that might be collapsed
        PHalfedge edge = collapse_candidate_queue.top();
        Halfedge_handle hh = edge.halfedge();
        collapse_candidate_queue.pop();
        collapse_candidate_queue.remove(PHalfedge(hh->opposite()));
        // step 2: try to collapse with the constraints of max_error
        Vertex_handle vh = collapse_applied(input_facet_tree,
            max_error_threshold, -1.0, true, NULL, &large_error_queue,
            &collapse_candidate_queue, hh, m_pRemesh);
        // step 3: update the adjacent halfedges
        if (vh != NULL) {
          ++nb_operations;
          if (m_track_information) {
            std::cout << "1 edge collapsed";
          }
        }
        if (m_track_information) {
          std::cout << std::endl;
        }
      }
    }
    std::cout << "Done (" << nb_operations << " local operations applied, "
              << timer.time() << " s)" << std::endl;
  }

  void split_local_longest_edge(const Facet_tree &input_facet_tree,
      const Bbox &bbox, Polyhedron *m_pRemesh) const {
    FT max_error = 0.0, min_radian = CGAL_PI;
    Halfedge_handle max_error_halfedge, min_radian_halfedge;
    FT max_error_threshold = get_max_error_threshold_value(bbox);
    max_error_halfedge = m_pRemesh->get_maximal_error(&max_error);
    min_radian_halfedge = m_pRemesh->get_minimal_radian(&min_radian);
    std::cout << "Split local longest edge..." << std::endl;
    std::cout << "(max error threshold = " << max_error_threshold
              << ", min angle threshold = " << m_min_angle_threshold
              << "бу)" << std::endl;
    Facet_handle fh = min_radian_halfedge->facet();
    Halfedge_handle longest_hh = m_pRemesh->get_longest_halfedge(fh);
    longest_hh = m_pRemesh->longest_side_propagation(longest_hh);
    Vertex_handle vh = split_edge(input_facet_tree, max_error_threshold,
        max_error, min_radian, false, NULL, NULL, longest_hh, m_pRemesh);
    if (vh != NULL) {
      std::cout << "1 local longest edge splitted" << std::endl;
    }
    else {
      std::cout << "Error: no edge splitted" << std::endl;
    }
  }

  void increase_minimal_angle(const Facet_tree &input_facet_tree,
                              const Bbox &bbox, Polyhedron *m_pRemesh) {
    // step 1: try to decrease the max error if necessary
    FT max_error_threshold = get_max_error_threshold_value(bbox);
    std::cout << "Increase minimal angle..." << std::endl;
    if (m_decrease_max_errors) {
      FT max_error = 0.0;
      Halfedge_handle max_error_halfedge;
      max_error_halfedge = m_pRemesh->get_maximal_error(&max_error);
      if (max_error >= max_error_threshold) {
        std::cout << "(max error threshold = " << max_error_threshold
          << ", max error = " << max_error << ")" << std::endl;
        greedy_reduce_error(input_facet_tree, max_error_threshold, max_error,
          true, false, NULL, NULL, max_error_halfedge,
          m_pRemesh);
        return;
      }
    }
    // step 2: try to increase the min radian
    FT min_radian = CGAL_PI;
    Halfedge_handle min_radian_halfedge;
    min_radian_halfedge = m_pRemesh->get_minimal_radian(&min_radian);
    std::cout << "(min angle threshold = " << m_min_angle_threshold
      << "бу, min angle = " << m_pRemesh->to_angle(min_radian)
      << "бу)" << std::endl;
    greedy_improve_angle(input_facet_tree, max_error_threshold, min_radian,
      true, NULL, NULL, min_radian_halfedge, m_pRemesh);
  }

  void maximize_minimal_angle(const Facet_tree &input_facet_tree,
                              const Bbox &bbox, Polyhedron *m_pRemesh) {
    FT max_error_threshold = get_max_error_threshold_value(bbox);
    FT min_radian_threshold = m_pRemesh->to_radian(m_min_angle_threshold);
    CGAL::Timer timer;
    timer.start();
    std::cout << "Greedy angles improvement..." << std::endl;
    std::cout << "(max error threshold = " << max_error_threshold
              << ", min angle threshold = " << m_min_angle_threshold
              << "бу, max mesh complexity = " << m_max_mesh_complexity
              << ")" << std::endl;
    FT max_error = 0, min_radian = CGAL_PI;
    Halfedge_handle max_error_halfedge, min_radian_halfedge;
    unsigned int nb_operations = 0;
    // Version 1: use dynamic priority queue
    DPQueue_halfedge_long large_error_queue;
    DPQueue_halfedge_short small_radian_queue;
    fill_small_radian_edges(max_error_threshold, &large_error_queue,
      &small_radian_queue, m_pRemesh);
    while (m_pRemesh->size_of_vertices() < m_max_mesh_complexity &&
      (!large_error_queue.empty() || !small_radian_queue.empty())) {
      while (m_pRemesh->size_of_vertices() < m_max_mesh_complexity &&
        !large_error_queue.empty()) {
        PHalfedge edge = large_error_queue.top();
        large_error_queue.pop();
        max_error_halfedge = edge.halfedge();
        max_error = CGAL::sqrt(edge.priority());
        if (m_track_information) {
          std::cout << ++nb_operations << ": max error = "
            << max_error << " ";
        }
        greedy_reduce_error(input_facet_tree, max_error_threshold,
          max_error, m_track_information, false, &large_error_queue,
          &small_radian_queue, max_error_halfedge, m_pRemesh);
      }
      if (m_pRemesh->size_of_vertices() < m_max_mesh_complexity &&
        !small_radian_queue.empty()) {
        PHalfedge edge = small_radian_queue.top();
        small_radian_queue.pop();
        min_radian_halfedge = edge.halfedge();
        min_radian = edge.priority();
        if (m_track_information) {
          std::cout << ++nb_operations << ": min angle = "
            << m_pRemesh->to_angle(min_radian) << " ";
        }
        greedy_improve_angle(input_facet_tree, max_error_threshold,
          min_radian, m_track_information, &large_error_queue,
          &small_radian_queue, min_radian_halfedge, m_pRemesh);
      }
    }
    // Version 2: do not use dynamic priority queue
    /*if (m_decrease_max_errors) {
      max_error_halfedge = m_pRemesh->get_maximal_error(&max_error);
      min_radian_halfedge = m_pRemesh->get_minimal_radian(&min_radian);
      while (m_pRemesh->size_of_vertices() < m_max_mesh_complexity &&
        (max_error >= max_error_threshold ||
        min_radian < min_radian_threshold)) {
        while (m_pRemesh->size_of_vertices() < m_max_mesh_complexity &&
          max_error >= max_error_threshold) {
          if (m_track_information) {
            std::cout << ++nb_operations << ": max error = "
              << max_error << " ";
          }
          greedy_reduce_error(input_facet_tree, max_error_threshold,
            max_error, m_track_information, false, NULL, NULL,
            max_error_halfedge, m_pRemesh);
          max_error_halfedge = m_pRemesh->get_maximal_error(&max_error);
          min_radian_halfedge = m_pRemesh->get_minimal_radian(&min_radian);
        }
        if (m_pRemesh->size_of_vertices() < m_max_mesh_complexity &&
          min_radian < min_radian_threshold) {
          if (m_track_information) {
            std::cout << ++nb_operations << ": min angle = "
              << m_pRemesh->to_angle(min_radian) << " ";
          }
          greedy_improve_angle(input_facet_tree, max_error_threshold,
            min_radian, m_track_information, NULL, NULL,
            min_radian_halfedge, m_pRemesh);
          max_error_halfedge = m_pRemesh->get_maximal_error(&max_error);
          min_radian_halfedge = m_pRemesh->get_minimal_radian(&min_radian);
        }
      }
    }
    else {
      min_radian_halfedge = m_pRemesh->get_minimal_radian(&min_radian);
      while (m_pRemesh->size_of_vertices() < m_max_mesh_complexity &&
        min_radian < min_radian_threshold) {
        if (m_track_information) {
          std::cout << ++nb_operations << ": min angle = "
            << m_pRemesh->to_angle(min_radian) << " ";
        }
        greedy_improve_angle(input_facet_tree, max_error_threshold,
          min_radian, m_track_information, NULL, NULL,
          min_radian_halfedge, m_pRemesh);
        min_radian_halfedge = m_pRemesh->get_minimal_radian(&min_radian);
      }
    }*/
    std::cout << "Done (" << nb_operations << " local operations applied, "
              << timer.time() << " s)" << std::endl;
  }

  void final_vertex_relocation(const Facet_tree &input_facet_tree,
                               const Bbox &bbox, Polyhedron *m_pRemesh) const {
    FT max_error_threshold = get_max_error_threshold_value(bbox);
    CGAL::Timer timer;
    timer.start();
    std::cout << "Final vertex relocation..." << std::endl;
    std::cout << "(max error threshold = " << max_error_threshold
      << ", min angle threshold = " << m_min_angle_threshold
      << "бу)..." << std::endl;
    DPQueue_vertex_short relocate_candidate_queue;
    fill_relocate_candidate_vertices(&relocate_candidate_queue, m_pRemesh);
    unsigned int index = 0, nb_relocate = 0;
    while (!relocate_candidate_queue.empty()) {
      if (m_track_information) {
        std::cout << ++index << ": relocate queue size = "
          << relocate_candidate_queue.size() << " ";
      }
      // step 1: get the top vertex that might be relocated
      PVertex vertex = relocate_candidate_queue.top();
      Vertex_handle vh = vertex.vertex();
      FT min_radian = vertex.priority();
      relocate_candidate_queue.pop();
      // step 2: try to relocate with the constrait of max_error and min_radian
      Point initial_point = get_initial_point_for_relocate(input_facet_tree,
        m_inherit_element_types, m_max_samples_per_area,
        m_feature_control_delta, vh, m_pRemesh);
      bool relocated = relocate_vertex(input_facet_tree, max_error_threshold,
        -1.0, min_radian, false, NULL, NULL, initial_point, vh, m_pRemesh);
      if (relocated) {
        update_relocate_candidate_vertices(vh, &relocate_candidate_queue,
          m_pRemesh);
        ++nb_relocate;
        if (m_track_information) {
          std::cout << "1 vertices relocated";
        }
      }
      if (m_track_information) {
        std::cout << std::endl;
      }
    }
    std::cout << "Done (" << nb_relocate << " vertices relocated, "
      << timer.time() << " s)" << std::endl;
  }

  // trees
  void build_facet_tree(const Polyhedron &polyhedron, const std::string &name,
                        Facet_tree *facet_tree) const {
    CGAL::Timer timer;
    timer.start();
    std::cout << "Build " << name << " facet AABB tree...";
    facet_tree->rebuild(faces(polyhedron).first,
      faces(polyhedron).second, polyhedron);
    facet_tree->accelerate_distance_queries();
    std::cout << "done (" << timer.time() << " s)" << std::endl;
  }

  // feature intensities
  void calculate_feature_intensities(const std::string &name,
                                     Polyhedron *polyhedron) const {
    polyhedron->calculate_feature_intensities("Input", m_dihedral_theta,
      m_dihedral_delta, m_sum_theta, m_sum_delta, m_inherit_element_types,
      m_feature_control_delta);
  }

  // links
  void generate_links(const Facet_tree &input_facet_tree,
      const Facet_tree &remesh_facet_tree, Polyhedron *m_pInput,
      Polyhedron *m_pRemesh) const {
    // step 1: clear all the links
    clear_links(m_pInput, m_pRemesh);
    // step 2: generate the out and in links
    generate_out_links(input_facet_tree, m_pRemesh);
    generate_in_links(remesh_facet_tree, m_pInput);
    // step 3: compute the max_squared_errors
    m_pRemesh->calculate_max_squared_errors();
  }

  void clear_links(Polyhedron *m_pInput, Polyhedron *m_pRemesh) const {
    // step 1: clear the out links
    m_pRemesh->clear_out_links();
    // step 2: clear the in links (out links from the perspective of m_pInput)
    m_pRemesh->clear_in_link_iterators();
    m_pInput->clear_out_links();
  }

  // visual elements
  void compute_facets(DrawType draw_type, RenderType render_type,
      const Bbox &bbox, bool is_input, Polyhedron *polyhedron,
      std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
      std::vector<float> *pos_face_colors, std::vector<float> *pos_boundaries,
      std::vector<float> *pos_samples) const {
    FT sum_theta_value = m_sum_theta * CGAL_PI;
    FT dihedral_theta_value = m_dihedral_theta * CGAL_PI;
    FT max_error_threshold_value = get_max_error_threshold_value(bbox);
    polyhedron->compute_facets(draw_type, render_type, m_inherit_element_types,
      m_feature_control_delta, sum_theta_value, dihedral_theta_value,
      max_error_threshold_value, is_input, pos_faces, pos_face_normals,
      pos_face_colors, pos_boundaries, pos_samples);
  }

 private:
  // max error
  FT get_max_error_threshold_value(const Bbox &bbox) const {
     FT diagonal = std::sqrt(std::pow(bbox.xmax() - bbox.xmin(), 2) +
       std::pow(bbox.ymax() - bbox.ymin(), 2) +
       std::pow(bbox.zmax() - bbox.zmin(), 2));
     return diagonal * m_max_error_threshold * 0.01;
  }

  // links
  void generate_out_links(const Facet_tree &input_facet_tree,
                          Polyhedron *m_pRemesh) const {
    // step 1: collect the edges and facets
    Facet_list facets;
    Halfedge_list edges;
    for (Facet_iterator fi = m_pRemesh->facets_begin();
      fi != m_pRemesh->facets_end(); ++fi) {
      facets.push_back(fi);
    }
    for (Edge_iterator ei = m_pRemesh->edges_begin();
      ei != m_pRemesh->edges_end(); ++ei) {
      edges.push_back(ei);
    }
    // step 2: calculate the number of samples per facet
    int nb_samples_per_facet_out = m_samples_per_facet_out;
    if (m_sample_number_strategy == SampleNumberStrategy::k_variable) {
      FT value = static_cast<double>(m_samples_per_facet_out);
      value *= input_facet_tree.size();
      value /= m_pRemesh->size_of_facets();
      nb_samples_per_facet_out = static_cast<int>(value);
    }
    m_pRemesh->get_nb_samples_per_facet(nb_samples_per_facet_out,
      m_max_samples_per_area,
      m_min_samples_per_triangle,
      m_sample_strategy,
      facets);
    // step 3: generating links
    CGAL::Timer timer;
    timer.start();
    std::cout << "Generating edge out links...";
    generate_edge_links(input_facet_tree, false, edges, m_pRemesh);
    std::cout << "Done, count: " << m_pRemesh->get_edge_out_link_count()
      << " (" << timer.time() << " s)" << std::endl;

    timer.reset();
    std::cout << "Generating vertex out links...";
    generate_vertex_links(input_facet_tree, false, m_pRemesh);
    std::cout << "Done, count: " << m_pRemesh->get_vertex_out_link_count()
      << " (" << timer.time() << " s)" << std::endl;

    timer.reset();
    std::cout << "Generating facet out links...";
    generate_facet_links(input_facet_tree, false, 5, facets, m_pRemesh);
    std::cout << "Done, count: " << m_pRemesh->get_facet_out_link_count()
      << " (" << timer.time() << " s)" << std::endl;

    m_pRemesh->reset_facet_tags(0, facets);
  }

  void generate_in_links(const Facet_tree &remesh_facet_tree,
                         Polyhedron *m_pInput) const {
    // step 1: collect the edges and facets
    Facet_list facets;
    Halfedge_list edges;
    for (Facet_iterator fi = m_pInput->facets_begin();
      fi != m_pInput->facets_end(); ++fi) {
      facets.push_back(fi);
    }
    for (Edge_iterator ei = m_pInput->edges_begin();
      ei != m_pInput->edges_end(); ++ei) {
      edges.push_back(ei);
    }
    // step 2: calculate the number of samples per facet
    m_pInput->get_nb_samples_per_facet(m_samples_per_facet_in,
      m_max_samples_per_area,
      m_min_samples_per_triangle,
      m_sample_strategy,
      facets);
    // step 3: generating links
    CGAL::Timer timer;
    timer.start();
    std::cout << "Generating edge in links...";
    generate_edge_links(remesh_facet_tree, true, edges, m_pInput);
    std::cout << "Done, count: " << m_pInput->get_edge_out_link_count()
      << " (" << timer.time() << " s)" << std::endl;

    timer.reset();
    std::cout << "Generating vertex in links...";
    generate_vertex_links(remesh_facet_tree, true, m_pInput);
    std::cout << "Done, count: " << m_pInput->get_vertex_out_link_count()
      << " (" << timer.time() << " s)" << std::endl;

    timer.reset();
    std::cout << "Generating facet in links...";
    generate_facet_links(remesh_facet_tree, true, 5, facets, m_pInput);
    std::cout << "Done, count: " << m_pInput->get_facet_out_link_count()
      << " (" << timer.time() << " s)" << std::endl;

    m_pInput->reset_facet_tags(0, facets);
  }

  void generate_local_links(const Facet_tree &input_facet_tree,
      bool reset_normal_dihedral, Link_iter_list *facet_in_links,
      Link_iter_list *edge_in_links, Link_pointer_list *vertex_in_links,
      Halfedge_handle hh, std::set<Facet_handle> *in_link_facets,
      Polyhedron *polyhedron) const {
    // step 1: update the local feature_intensity around vh
    Vertex_handle vp = polyhedron->get_source_vertex(hh);
    Vertex_handle vq = polyhedron->get_target_vertex(hh);
    polyhedron->update_local_feature_intensity(vp, reset_normal_dihedral,
      m_dihedral_theta, m_dihedral_delta, m_sum_theta, m_sum_delta);
    polyhedron->update_local_feature_intensity(vp, reset_normal_dihedral,
      m_dihedral_theta, m_dihedral_delta, m_sum_theta, m_sum_delta);
    // step 2: clear local links
    polyhedron->clear_local_links(hh, in_link_facets);
    // step 3: generate local out links
    generate_local_out_links(input_facet_tree, hh, polyhedron);
    // step 4: generate local in links
    generate_local_in_links(in_link_facets,
      facet_in_links, edge_in_links, vertex_in_links, polyhedron);
  }

  void generate_local_links(const Facet_tree &input_facet_tree,
      bool reset_normal_dihedral, Link_iter_list *facet_in_links,
      Link_iter_list *edge_in_links, Link_pointer_list *vertex_in_links,
      Vertex_handle vh, std::set<Facet_handle> *in_link_facets,
      Polyhedron *polyhedron) const {
    // step 1: update the local feature_intensity aournd vh
    polyhedron->update_local_feature_intensity(vh, reset_normal_dihedral,
      m_dihedral_theta, m_dihedral_delta, m_sum_theta, m_sum_delta);
    // step 2: clear local links
    polyhedron->clear_local_links(vh, in_link_facets);
    // step 3: generate local out links
    generate_local_out_links(input_facet_tree, vh, polyhedron);
    // step 4: generate local in links
    generate_local_in_links(in_link_facets,
      facet_in_links, edge_in_links, vertex_in_links, polyhedron);
  }

  void generate_local_out_links(const Facet_tree &input_facet_tree,
      Vertex_handle vh, Polyhedron *polyhedron) const {
    // step 1: collect the one_ring facets and edges
    Facet_list facets;
    Halfedge_list edges;
    Halfedge_around_vertex_circulator vcirc = vh->vertex_begin();
    Halfedge_around_vertex_circulator vend = vcirc;
    CGAL_For_all(vcirc, vend) {
      if (!vcirc->is_border()) {
        facets.push_back(vcirc->facet());
      }
      edges.push_back(vcirc);
    }
    // step 2: calculate the number of samples per facets
    polyhedron->get_nb_samples_per_facet(m_samples_per_facet_out,
      m_max_samples_per_area,
      m_min_samples_per_triangle,
      m_sample_strategy,
      facets);
    // step 3: generate out links (is_in_link is set to false)
    generate_edge_links(input_facet_tree, false, edges, polyhedron);
    generate_vertex_link(input_facet_tree, false, vh, polyhedron);
    generate_facet_links(input_facet_tree, false, m_bvd_iteration_count,
      facets, polyhedron);
    // step 4: reset the facet tags
    polyhedron->reset_facet_tags(0, facets);
  }

  void generate_local_out_links(const Facet_tree &input_facet_tree,
      Halfedge_handle hh, Polyhedron *polyhedron) const {
    // step 1: collect the facets and edges
    Facet_list facets;
    Halfedge_list edges;
    if (!hh->is_border()) {
      facets.push_back(hh->facet());
    }
    if (!hh->opposite()->is_border()) {
      facets.push_back(hh->opposite()->facet());
    }
    edges.push_back(hh);
    // step 2: calculate the number of samplers per facets
    polyhedron->get_nb_samples_per_facet(m_samples_per_facet_out,
      m_max_samples_per_area,
      m_min_samples_per_triangle,
      m_sample_strategy,
      facets);
    // step 3: generate out links (is_in_link is set to false)
    generate_edge_links(input_facet_tree, false, edges, polyhedron);
    generate_facet_links(input_facet_tree, false, m_bvd_iteration_count,
      facets, polyhedron);
    // step 4: reset the facet tags
    polyhedron->reset_facet_tags(0, facets);
  }

  void generate_local_in_links(std::set<Facet_handle> *in_link_facets,
      Link_iter_list *facet_in_links, Link_iter_list *edge_in_links,
      Link_pointer_list *vertex_in_links, Polyhedron *polyhedron) const {
    Facet_tree::Point_and_primitive_id pp;
    Facet_handle closest_facet;
    if (m_use_local_aabb_tree) {
      // step 1: build the local aabb tree
      Facet_tree local_facet_tree;
      local_facet_tree.rebuild(
        in_link_facets->begin(), in_link_facets->end(), *polyhedron);
      local_facet_tree.accelerate_distance_queries();
      // step 2: update the in links
      for (Link_iter_list_iter it = facet_in_links->begin();
        it != facet_in_links->end(); ++it) {
        Link_list_iter llit = *it;
        pp = local_facet_tree.closest_point_and_primitive(llit->second.first);
        llit->second.second = pp.first;   // update the closest point
        closest_facet = pp.second;
        closest_facet->facet_in_links().push_back(llit);
      }
      for (Link_iter_list_iter it = edge_in_links->begin();
        it != edge_in_links->end(); ++it) {
        Link_list_iter llit = *it;
        pp = local_facet_tree.closest_point_and_primitive(llit->second.first);
        llit->second.second = pp.first;   // update the closest point
        closest_facet = pp.second;
        closest_facet->edge_in_links().push_back(llit);
      }
      for (Link_pointer_iter it = vertex_in_links->begin();
        it != vertex_in_links->end(); ++it) {
        Link *link = *it;
        pp = local_facet_tree.closest_point_and_primitive(link->second.first);
        link->second.second = pp.first;   // update the closest point
        closest_facet = pp.second;
        closest_facet->vertex_in_links().push_back(link);
      }
    }
    else {
      // for each sample, update its closest point, push back to the primitive
      for (Link_iter_list_iter it = facet_in_links->begin();
        it != facet_in_links->end(); ++it) {
        Link_list_iter llit = *it;
        pp = get_closest_point_and_primitive(*in_link_facets,
          llit->second.first, polyhedron);
        llit->second.second = pp.first;   // update the cloest point
        closest_facet = pp.second;
        closest_facet->facet_in_links().push_back(llit);
      }
      for (Link_iter_list_iter it = edge_in_links->begin();
        it != edge_in_links->end(); ++it) {
        Link_list_iter llit = *it;
        pp = get_closest_point_and_primitive(*in_link_facets,
          llit->second.first, polyhedron);
        llit->second.second = pp.first;   // update the closest point
        closest_facet = pp.second;
        closest_facet->edge_in_links().push_back(llit);
      }
      for (Link_pointer_iter it = vertex_in_links->begin();
        it != vertex_in_links->end(); ++it) {
        Link *link = *it;
        pp = get_closest_point_and_primitive(*in_link_facets,
          link->second.first, polyhedron);
        link->second.second = pp.first;   // update the closest point
        closest_facet = pp.second;
        closest_facet->vertex_in_links().push_back(link);
      }
    }
  }

  void generate_edge_links(const Facet_tree &facet_tree,
      bool is_in_link, const Halfedge_list &edges,
      Polyhedron *polyhedron) const {
    // precondition: the number of samples per facets has been calculated
    for (auto it = edges.begin(); it != edges.end(); ++it) {
      Halfedge_handle hh = *it;
      if (hh->normal_dihedral() == -1.0) {  // value stored in its opposite
        hh = hh->opposite();
      }
      FT facet_area = polyhedron->area(hh->facet());
      int nb_edge_out_links = hh->facet()->tag();
      if (!hh->opposite()->is_border()) {
        facet_area += polyhedron->area(hh->opposite()->facet());
        nb_edge_out_links += hh->opposite()->facet()->tag();
      }
      if (nb_edge_out_links > 0) {
        FT area_per_sample = facet_area / nb_edge_out_links;
        FT diameter = 2 * CGAL::sqrt(area_per_sample / CGAL_PI);
        int nb_samples = polyhedron->length(hh) / diameter;
        nb_samples = std::max(nb_samples, m_min_samples_per_triangle);
        FT capacity =
          m_use_stratified_sampling ? facet_area / (3 * nb_samples) : 1.0;
        sample_edge_links(facet_tree, is_in_link, capacity, nb_samples, hh);
      }
    }
  }

  void sample_edge_links(const Facet_tree &facet_tree, bool is_in_link,
      FT capacity, int nb_samples, Halfedge_handle hh) const {
    Vertex_handle vp = hh->opposite()->vertex(), vq = hh->vertex();
    Vector vector = vq->point() - vp->point();
    FT step = 1.0 / (nb_samples + 1);
    Facet_tree::Point_and_primitive_id pp;
    for (int i = 1; i <= nb_samples; ++i) {
      const Point sample = vp->point() + vector * step * i;
      pp = facet_tree.closest_point_and_primitive(sample);
      FT feature_weight = vp->feature_intensity() * (nb_samples + 1 - i) +
        vq->feature_intensity() * i;
      feature_weight /= (nb_samples + 1);   // interpolated feature intensity
      // 1) insert the sample in the source
      Link_list_iter it =
        hh->edge_out_links().insert(hh->edge_out_links().end(),
        std::make_pair(feature_weight * capacity,
        std::make_pair(sample, pp.first)));
      // 2) insert the iterator in the target if necessary
      if (is_in_link) {
        Facet_handle closest_facet = pp.second;
        closest_facet->edge_in_links().push_back(it);
      }
    }
  }

  void generate_vertex_links(const Facet_tree &facet_tree, bool is_in_link,
                             Polyhedron *polyhedron) const {
    for (Vertex_iterator vi = polyhedron->vertices_begin();
      vi != polyhedron->vertices_end(); ++vi) {
      generate_vertex_link(facet_tree, is_in_link, vi, polyhedron);
    }
  }

  void generate_vertex_link(const Facet_tree &facet_tree, bool is_in_link,
                            Vertex_handle vh, Polyhedron *polyhedron) const {
    // precondition: the weights of the vertices have been calcualted
    Facet_tree::Point_and_primitive_id pp;
    FT capacity =
      m_use_stratified_sampling ? polyhedron->get_vertex_capacity(vh) : 1.0;
    pp = facet_tree.closest_point_and_primitive(vh->point());
    // 1) insert the sample in the source
    vh->vertex_out_link().first = capacity * vh->feature_intensity();
    vh->vertex_out_link().second.first = vh->point();
    vh->vertex_out_link().second.second = pp.first;
    // 2) insert the sample in the target if necessary
    if (is_in_link) {
      Facet_handle closest_facet = pp.second;
      closest_facet->vertex_in_links().push_back(&vh->vertex_out_link());
    }
  }

  void generate_facet_links(const Facet_tree &facet_tree, bool is_in_link,
      int bvd_iteration_count, const Facet_list &facets,
      Polyhedron *polyhedron) const {
    Point_list inner_samples;
    Point_iter pit;
    std::list<double> feature_weights;
    std::list<double>::iterator fit;
    Facet_tree::Point_and_primitive_id pp;
    for (auto it = facets.begin(); it != facets.end(); ++it) {
      Facet_handle fh = *it;
      if (fh->tag() > 0) {
        polyhedron->generate_random_samples(m_use_stratified_sampling,
          fh->tag(), bvd_iteration_count, fh,
          &inner_samples, &feature_weights);
        FT capacity =
          m_use_stratified_sampling ? polyhedron->area(fh) / fh->tag() : 1.0;
        for (pit = inner_samples.begin(), fit = feature_weights.begin();
          pit != inner_samples.end(); ++pit, ++fit) {
          pp = facet_tree.closest_point_and_primitive(*pit);
          // 1) insert the sample in the source
          Link_list_iter it =
            fh->facet_out_links().insert(fh->facet_out_links().end(),
            std::make_pair(capacity * (*fit),
            std::make_pair(*pit, pp.first)));
          // 2) insert the samplein the target if necessary
          if (is_in_link) {
            Facet_handle closest_facet = pp.second;
            closest_facet->facet_in_links().push_back(it);
          }
        }
      }
    }
  }

  void backup_local_links(const std::set<Facet_handle> &extended_facets,
      std::map<Facet_handle, Link_list> *facet_out_map,
      std::map<Halfedge_handle, Link_list> *edge_out_map,
      Link *vertex_out, Link_iter_list *facet_in_links,
      Link_iter_list *edge_in_links, Link_pointer_list *vertex_in_links,
      Vertex_handle vh) const {
    // step 1: backup the out links
    Halfedge_around_vertex_circulator vcirc = vh->vertex_begin();
    Halfedge_around_vertex_circulator vend = vcirc;
    CGAL_For_all(vcirc, vend) {
      if (!vcirc->is_border()) {
        Facet_handle fh = vcirc->facet();
        (*facet_out_map)[fh].insert((*facet_out_map)[fh].end(),
          fh->facet_out_links().begin(), fh->facet_out_links().end());
      }
      Halfedge_handle hh = vcirc;
      if (hh->normal_dihedral() == -1.0) {
        hh = hh->opposite();
      }
      (*edge_out_map)[hh].insert((*edge_out_map)[hh].end(),
        hh->edge_out_links().begin(), hh->edge_out_links().end());
    }
    (*vertex_out) = vh->vertex_out_link();
    // step 2: backup the in links
    backup_local_in_links(extended_facets, facet_in_links,
      edge_in_links, vertex_in_links);
  }

  void backup_local_in_links(const std::set<Facet_handle> &extended_facets,
      Link_iter_list *facet_in_links, Point_list *facet_in_end_points,
      Link_iter_list *edge_in_links, Point_list *edge_in_end_points,
      Link_pointer_list *vertex_in_links,
      Point_list *vertex_in_end_points) const {
    // step 1: store the in link iterations
    backup_local_in_links(extended_facets, facet_in_links,
      edge_in_links, vertex_in_links);
    // step 2: store the in lin end points (for simulation purpose)
    for (Link_iter_list_iter it = facet_in_links->begin();
      it != facet_in_links->end(); ++it) {
      Link_list_iter llit = *it;
      facet_in_end_points->push_back(llit->second.second);
    }
    for (Link_iter_list_iter it = edge_in_links->begin();
      it != edge_in_links->end(); ++it) {
      Link_list_iter llit = *it;
      edge_in_end_points->push_back(llit->second.second);
    }
    for (Link_pointer_iter it = vertex_in_links->begin();
      it != vertex_in_links->end(); ++it) {
      Link* link = *it;
      vertex_in_end_points->push_back(link->second.second);
    }
  }

  void backup_local_in_links(const std::set<Facet_handle> &extended_facets,
      Link_iter_list *facet_in_links, Link_iter_list *edge_in_links,
      Link_pointer_list *vertex_in_links) const {
    // step 2: backup in links
    for (typename std::set<Facet_handle>::iterator it = extended_facets.begin();
      it != extended_facets.end(); ++it) {
      Facet_const_handle fh = *it;
      // store facet in links
      facet_in_links->insert(facet_in_links->end(),
        fh->facet_in_links().begin(), fh->facet_in_links().end());
      // store edge in links
      edge_in_links->insert(edge_in_links->end(),
        fh->edge_in_links().begin(), fh->edge_in_links().end());
      // store vertex in links
      vertex_in_links->insert(vertex_in_links->end(),
        fh->vertex_in_links().begin(), fh->vertex_in_links().end());
    }
  }

  void restore_local_links(const Point_list &facet_in_end_points,
      Link_iter_list *facet_in_links, const Point_list &edge_in_end_points,
      Link_iter_list *edge_in_links, const Point_list &vertex_in_end_points,
      Link_pointer_list *vertex_in_links) const {
    Link_iter_list_iter lit;
    Point_const_iter pit;
    Link_pointer_iter link_iter;
    for (lit = facet_in_links->begin(), pit = facet_in_end_points.begin();
      lit != facet_in_links->end(); ++lit, ++pit) {
      Link_list_iter it = *lit;
      it->second.second = *pit;
    }
    for (lit = edge_in_links->begin(), pit = edge_in_end_points.begin();
      lit != edge_in_links->end(); ++lit, ++pit) {
      Link_list_iter it = *lit;
      it->second.second = *pit;
    }
    for (link_iter = vertex_in_links->begin(),
      pit = vertex_in_end_points.begin(); link_iter != vertex_in_links->end();
      ++link_iter, ++pit) {
      Link *link = *link_iter;
      link->second.second = *pit;
    }
  }

  void restore_local_links(
      const std::map<Facet_handle, Link_list> &facet_out_map,
      const std::map<Halfedge_handle, Link_list> &edge_out_map,
      const Link &vertex_out, Link_iter_list *facet_in_links,
      Link_iter_list *edge_in_links, Link_pointer_list *vertex_in_links,
      Vertex_handle vh, std::set<Facet_handle> *in_link_facets,
      Polyhedron *polyhedron) const {
    // step 1: update the local feature_intensity around vh
    polyhedron->update_local_feature_intensity(vh, false,
      m_dihedral_theta, m_dihedral_delta, m_sum_theta, m_sum_delta);
    // step 2: clear local_links
    polyhedron->clear_local_links(vh, in_link_facets);
    // step 3: restore local out links
    restore_local_out_links(facet_out_map, edge_out_map, vertex_out, vh);
    // step 4: generate local in links
    generate_local_in_links(in_link_facets,
      facet_in_links, edge_in_links, vertex_in_links, polyhedron);
  }

  void restore_local_out_links(
      const std::map<Facet_handle, Link_list> &facet_out_map,
      const std::map<Halfedge_handle, Link_list> &edge_out_map,
      const Link &vertex_out, Vertex_handle vh) const {
    // step 1: restore the facet out links
    for (auto it = facet_out_map.begin(); it != facet_out_map.end(); ++it) {
      Facet_handle fh = it->first;
      fh->facet_out_links().insert(fh->facet_out_links().end(),
        it->second.begin(), it->second.end());
    }
    // step 2: restore the edge out links
    for (auto it = edge_out_map.begin(); it != edge_out_map.end(); ++it) {
      Halfedge_handle hh = it->first;
      hh->edge_out_links().insert(hh->edge_out_links().end(),
        it->second.begin(), it->second.end());
    }
    vh->vertex_out_link() = vertex_out;
  }

  // queues (for eliminating degenerated facets)
  void fill_degenerated_facets_queue(FT radian_threshold,
      DPQueue_facet_long *queue, Polyhedron *polyhedron) const {
    for (Facet_iterator fi = polyhedron->facets_begin();
      fi != polyhedron->facets_end(); ++fi) {
      if (polyhedron->get_smallest_radian(fi) < radian_threshold) {
        Halfedge_handle hh = polyhedron->get_longest_halfedge(fi);
        queue->push(PFacet(fi, polyhedron->squared_length(hh)));
      }
    }
  }

  void remove_incident_facets(Halfedge_handle hh,
      DPQueue_facet_long *queue, Polyhedron *polyhedron) const {
    Vertex_handle vp = polyhedron->get_source_vertex(hh);
    Vertex_handle vq = polyhedron->get_target_vertex(hh);
    Halfedge_around_vertex_circulator vcirc = vp->vertex_begin();
    Halfedge_around_vertex_circulator vend = vcirc;
    CGAL_For_all(vcirc, vend) {
      if (!vcirc->is_border()) {
        queue->remove(PFacet(vcirc->facet()));
      }
    }
    vcirc = vq->vertex_begin();
    vend = vcirc;
    CGAL_For_all(vcirc, vend) {
      if (!vcirc->is_border()) {
        queue->remove(PFacet(vcirc->facet()));
      }
    }
  }

  void add_circulator_degenerate_facets(Vertex_handle vh, FT radian_threshold,
      DPQueue_facet_long *queue, Polyhedron *polyhedron) const {
    Halfedge_around_vertex_circulator he = vh->vertex_begin();
    Halfedge_around_vertex_circulator end = he;
    FT longest_squared_length = 0.0;
    CGAL_For_all(he, end) {
      if (!he->is_border()) {
        Facet_handle fh = he->facet();
        Halfedge_handle shortest_hh = polyhedron->get_shortest_halfedge(fh);
        FT smallest_radian = polyhedron->get_opposite_radian(shortest_hh);
        if (smallest_radian < radian_threshold) {
          // degenerate edge
          if (polyhedron->squared_length(shortest_hh) < SQUARED_MIN_VALUE) {
            queue->push(PFacet(fh, MAX_VALUE));
          }
          // degenerate facet
          else {
            Halfedge_const_handle hh = polyhedron->get_longest_halfedge(fh);
            longest_squared_length = polyhedron->squared_length(hh);
            queue->push(PFacet(fh, longest_squared_length));
          }
        }
      }
    }
  }

  // queues (for splitting long edges)
  void fill_queue_with_long_edges(FT max_squared_length,
      DPQueue_halfedge_long *queue, Polyhedron *polyhedron) const {
    for (Edge_iterator ei = polyhedron->edges_begin();
      ei != polyhedron->edges_end(); ++ei) {
      const FT sl = polyhedron->squared_length(ei);
      if (sl >= max_squared_length) {
        Halfedge_handle hh = ei->is_border() ? ei->opposite() : ei;
        queue->push(PHalfedge(hh, sl));
      }
    }
  }

  void add_circular_long_edges(FT max_squared_length, Vertex_handle vh,
      DPQueue_halfedge_long *queue, Polyhedron *polyhedron) const {
    Halfedge_around_vertex_circulator vcirc = vh->vertex_begin();
    Halfedge_around_vertex_circulator vend = vcirc;
    CGAL_For_all(vcirc, vend) {
      const FT sl = polyhedron->squared_length(vcirc);
      if (sl >= max_squared_length) {
        Halfedge_handle hh = vcirc->is_border() ? vcirc->opposite() : vcirc;
        queue->push(PHalfedge(hh, sl));
      }
    }
  }

  // queues (for isotropic remeshing)
  void fill_collapse_candidate_edges(FT max_error_threshold,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *collapse_candidate_queue,
      Polyhedron *m_pRemesh) const {
    // step 1: fill the large error queue if necessary
    if (m_decrease_max_errors) {
      FT max_se_threshold = max_error_threshold * max_error_threshold;
      Halfedge_handle max_error_halfedge;
      for (Facet_iterator fi = m_pRemesh->facets_begin();
        fi != m_pRemesh->facets_end(); ++fi) {
        FT max_se = fi->max_squared_error();
        if (max_se >= max_se_threshold) {
          max_error_halfedge = m_pRemesh->get_longest_halfedge(fi);
          large_error_queue->push(PHalfedge(max_error_halfedge, max_se));
        }
      }
    }
    // step 2: fill the collapse_candidate queue
    for (Halfedge_iterator hi = m_pRemesh->halfedges_begin();
      hi != m_pRemesh->halfedges_end(); ++hi) {
      if (!hi->is_border()) {
        const FT length = m_pRemesh->length(hi);
        FT radian = m_pRemesh->get_opposite_radian(hi);
        collapse_candidate_queue->push(PHalfedge(hi, length * radian));
      }
    }
  }

  void fill_small_radian_edges(FT max_error_threshold,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_radian_queue,
      Polyhedron *m_pRemesh) const {
    // step 1: fill the large error queue if necessary
    if (m_decrease_max_errors) {
      FT max_se_threshold = max_error_threshold * max_error_threshold;
      Halfedge_handle max_error_halfedge;
      for (Facet_iterator fi = m_pRemesh->facets_begin();
        fi != m_pRemesh->facets_end(); ++fi) {
        FT max_se = fi->max_squared_error();
        if (max_se >= max_se_threshold) {
          max_error_halfedge = m_pRemesh->get_longest_halfedge(fi);
          large_error_queue->push(PHalfedge(max_error_halfedge, max_se));
        }
      }
    }
    // step 2: fill the small radian queue
    FT min_radian_threshold = m_pRemesh->to_radian(m_min_angle_threshold);
    for (Halfedge_iterator hi = m_pRemesh->halfedges_begin();
      hi != m_pRemesh->halfedges_end(); ++hi) {
      if (!hi->is_border()) {
        FT radian = m_pRemesh->get_opposite_radian(hi);
        if (radian < min_radian_threshold) {
          small_radian_queue->push(PHalfedge(hi, radian));
        }
      }
    }
  }

  void fill_relocate_candidate_vertices(
      DPQueue_vertex_short *relocate_candidate_queue,
      Polyhedron *m_pRemesh) const {
    for (Vertex_iterator vi = m_pRemesh->vertices_begin();
      vi != m_pRemesh->vertices_end(); ++vi) {
      FT min_radian = m_pRemesh->get_minimal_radian_incident_to_vertex(vi);
      relocate_candidate_queue->push(PVertex(vi, min_radian));
    }
  }

  void update_relocate_candidate_vertices(Vertex_handle vh,
      DPQueue_vertex_short *relocate_candidate_queue,
      Polyhedron *m_pRemesh) const {
    std::set<Vertex_handle> incident_vertices;
    m_pRemesh->collect_indicent_vertices(vh, &incident_vertices);
    incident_vertices.insert(vh);
    for (auto it = incident_vertices.begin();
      it != incident_vertices.end(); ++it) {
      Vertex_handle vh = *it;
      FT min_radian = m_pRemesh->get_minimal_radian_incident_to_vertex(vh);
      relocate_candidate_queue->push_or_update(PVertex(vh, min_radian));
    }
  }

  void remove_small_value_edges_before_collapse(Halfedge_handle hh,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue, Polyhedron *m_pRemesh) const {
    std::set<Facet_handle> one_ring_facets;
    m_pRemesh->collect_one_ring_facets_incident_to_edge(hh, &one_ring_facets);
    remove_small_value_edges(one_ring_facets, large_error_queue,
      small_value_queue, m_pRemesh);
  }

  void add_small_value_edges_after_collapse(Vertex_handle vh,
      FT max_error_threshold, bool reduce_complexity,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue, Polyhedron *m_pRemesh) const {
    std::set<Facet_handle> one_ring_facets;
    m_pRemesh->collect_one_ring_facets_incident_to_vertex(vh,
      &one_ring_facets);
    add_small_value_edges(one_ring_facets, max_error_threshold,
      reduce_complexity, large_error_queue, small_value_queue, m_pRemesh);
  }

  void remove_small_value_edges_before_split(Halfedge_handle hh,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue, Polyhedron *m_pRemesh) const {
    std::set<Facet_handle> one_ring_facets;
    m_pRemesh->collect_facets_incident_to_edge(hh, &one_ring_facets);
    remove_small_value_edges(one_ring_facets, large_error_queue,
      small_value_queue, m_pRemesh);
  }

  void add_small_value_edges_after_split(Vertex_handle vh,
      FT max_error_threshold, bool reduce_complexity,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue, Polyhedron *m_pRemesh) const {
    std::set<Facet_handle> one_ring_facets;
    m_pRemesh->collect_one_ring_facets_incident_to_vertex(vh,
      &one_ring_facets);
    add_small_value_edges(one_ring_facets, max_error_threshold,
      reduce_complexity, large_error_queue, small_value_queue, m_pRemesh);
  }

  void remove_small_value_edges_before_flip(Halfedge_handle hh,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue, Polyhedron *m_pRemesh) const {
    std::set<Facet_handle> one_ring_facets;
    m_pRemesh->collect_facets_incident_to_edge(hh, &one_ring_facets);
    remove_small_value_edges(one_ring_facets, large_error_queue,
      small_value_queue, m_pRemesh);
  }

  void add_small_value_edges_after_flip(Halfedge_handle hh,
      FT max_error_threshold, bool reduce_complexity,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue, Polyhedron *m_pRemesh) const {
    std::set<Facet_handle> one_ring_facets;
    m_pRemesh->collect_facets_incident_to_edge(hh, &one_ring_facets);
    add_small_value_edges(one_ring_facets, max_error_threshold,
      reduce_complexity, large_error_queue, small_value_queue, m_pRemesh);
  }

  void remove_small_value_edges_before_relocate(Vertex_handle vh,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue, Polyhedron *m_pRemesh) const {
    std::set<Facet_handle> one_ring_facets;
    m_pRemesh->collect_one_ring_facets_incident_to_vertex(vh,
      &one_ring_facets);
    remove_small_value_edges(one_ring_facets, large_error_queue,
      small_value_queue, m_pRemesh);
  }

  void add_small_value_edges_after_relocate(Vertex_handle vh,
      FT max_error_threshold, bool reduce_complexity,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue, Polyhedron *m_pRemesh) const {
    std::set<Facet_handle> one_ring_facets;
    m_pRemesh->collect_one_ring_facets_incident_to_vertex(vh,
      &one_ring_facets);
    add_small_value_edges(one_ring_facets, max_error_threshold,
      reduce_complexity, large_error_queue, small_value_queue, m_pRemesh);
  }

  void remove_small_value_edges(const std::set<Facet_handle> &one_ring_facets,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue, Polyhedron *m_pRemesh) const {
    // small_value_queue may be small_radian_queue or collapse_candidate_queue
    // step 1: remove from small_value_queue
    typename std::set<Facet_handle>::iterator it;
    for (it = one_ring_facets.begin(); it != one_ring_facets.end(); ++it) {
      Facet_handle fh = *it;
      Halfedge_handle hh = fh->halfedge();
      small_value_queue->remove(PHalfedge(hh));
      small_value_queue->remove(PHalfedge(hh->next()));
      small_value_queue->remove(PHalfedge(hh->prev()));
    }
    // step 2: remove from large_error_queue if necessary
    if (m_decrease_max_errors) {
      std::set<Facet_handle> extended_facets;
      m_pRemesh->extend_facets(one_ring_facets, m_stencil_ring_size,
        &extended_facets);
      for (it = extended_facets.begin(); it != extended_facets.end(); ++it) {
        Facet_handle fh = *it;
        Halfedge_handle hh = fh->halfedge();
        large_error_queue->remove(PHalfedge(hh));
        large_error_queue->remove(PHalfedge(hh->next()));
        large_error_queue->remove(PHalfedge(hh->prev()));
      }
    }
  }

  void add_small_value_edges(const std::set<Facet_handle> &one_ring_facets,
      FT max_error_threshold, bool reduce_complexity,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue, Polyhedron *m_pRemesh) const {
    // small_value_queue may be small_radian_queue or collapse_candidate_queue
    // If is former, the prority is radian; otherwise it is length * radian
    // step 1: add to small_value_queue
    typename std::set<Facet_handle>::iterator it;
    if (reduce_complexity) {    // for min rdian improvement
      for (it = one_ring_facets.begin(); it != one_ring_facets.end(); ++it) {
        Facet_handle fh = *it;
        Halfedge_handle hh = fh->halfedge();
        for (int i = 0; i <= 2; ++i) {
          FT radian = m_pRemesh->get_opposite_radian(hh);
          FT length = m_pRemesh->length(hh);
          small_value_queue->push(PHalfedge(hh, radian * length));
          hh = hh->next();
        }
      }
    }
    else {                       // for mesh complexity reduction
      FT min_radian_threshold = m_pRemesh->to_radian(m_min_angle_threshold);
      for (it = one_ring_facets.begin(); it != one_ring_facets.end(); ++it) {
        Facet_handle fh = *it;
        Halfedge_handle hh = fh->halfedge();
        for (int i = 0; i <= 2; ++i) {
          FT radian = m_pRemesh->get_opposite_radian(hh);
          if (radian < min_radian_threshold) {
            small_value_queue->push(PHalfedge(hh, radian));
          }
          hh = hh->next();
        }
      }
    }
    // step 2: add to large_error_queue if necessary
    if (m_decrease_max_errors) {
      FT max_se_threshold = max_error_threshold * max_error_threshold;
      std::set<Facet_handle> extended_facets;
      m_pRemesh->extend_facets(one_ring_facets, m_stencil_ring_size,
        &extended_facets);
      for (it = extended_facets.begin(); it != extended_facets.end(); ++it) {
        Facet_handle fh = *it;
        FT max_se = fh->max_squared_error();
        if (max_se >= max_se_threshold) {
          Halfedge_handle longest_hh = m_pRemesh->get_longest_halfedge(fh);
          large_error_queue->push(PHalfedge(longest_hh, max_se));
        }
      }
    }
  }

  // isotropic remeshing
  void greedy_improve_angle(const Facet_tree &input_facet_tree,
      FT max_error_threshold, FT min_radian, bool track_information,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_radian_queue,
      Halfedge_handle min_radian_halfedge, Polyhedron *m_pRemesh) {
    // improve min_radian constrained by max_error_threshold in the following:
    // 1) If collapse applies, collapse and return;
    // 2) If flip applies, flip and return;
    // 3) If relocate applies, relocate and return;
    // 4) Find the local_longest_hh. If flip applies, flip local_longest_hh;
    //                               Otherwwiese, split local_longest_hh.
    // step 1: try to collapse
    bool infinite_loop = false;
    Vertex_handle vh = collapse_applied(input_facet_tree, max_error_threshold,
      min_radian, false, &infinite_loop, large_error_queue,
      small_radian_queue, min_radian_halfedge, m_pRemesh);
    if (vh != NULL) {
      if (track_information) {
        std::cout << "1 edge collapsed" << std::endl;
      }
      return;
    }
    // if no infinite loop encountered, try flip and relocate
    if (!infinite_loop) {
      // step 2: try to flip
      if (m_apply_edge_flip) {
        int nb_flip = flip_applied(input_facet_tree, max_error_threshold,
          -1.0, min_radian, false, large_error_queue, small_radian_queue,
          min_radian_halfedge, m_pRemesh);
        if (nb_flip > 0) {
          if (track_information) {
            std::cout << nb_flip << " edges flipped" << std::endl;
          }
          return;
        }
      }
      // step 3: try to relocate
      int nb_relocate = relocate_applied(input_facet_tree,
        max_error_threshold, -1.0, min_radian, false, large_error_queue,
        small_radian_queue, min_radian_halfedge, m_pRemesh);
      if (nb_relocate > 0) {
        if (track_information) {
          std::cout << nb_relocate << " vertices relocated" << std::endl;
        }
        return;
      }
    }
    // step 3: split for later improvement
    Facet_handle fh = min_radian_halfedge->facet();
    Halfedge_handle longest_hh = m_pRemesh->get_longest_halfedge(fh);
    longest_hh = m_pRemesh->longest_side_propagation(longest_hh);
    if (m_apply_edge_flip) {
      Halfedge_handle hnew = flip_edge(input_facet_tree, max_error_threshold,
        -1.0, min_radian, false, large_error_queue, small_radian_queue,
        longest_hh, m_pRemesh);
      if (hnew != NULL) {
        if (track_information) {
          std::cout << "1 longest edge flipped" << std::endl;
        }
        return;
      }
    }
    split_edge(input_facet_tree, max_error_threshold, -1.0, min_radian,
      false, large_error_queue, small_radian_queue, longest_hh, m_pRemesh);
    if (track_information) {
      std::cout << "1 edge splitted" << std::endl;
    }
  }

  void greedy_reduce_error(const Facet_tree &input_facet_tree,
      FT max_error_threshold, FT max_error, bool track_information,
      bool reduce_complexity, DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue,
      Halfedge_handle max_error_halfedge, Polyhedron *m_pRemesh) const {
    // 1. reduce max_error in the following order:
    // 1) If flip reduces the max_error, flip and return;
    // 2) If relocate reduces the max_error, relocate and return;
    // 3) Otherwise, split the max_error_halfedge.
    // 2. if improve_min_radian is true, small_value_queue is the
    //    small_radian_queue; otherwise, it is the collapse_candidate_queue.
    // step 1: try to flip
    if (m_apply_edge_flip && flip_edge(input_facet_tree, -1.0, max_error,
      -1.0, reduce_complexity, large_error_queue, small_value_queue,
      max_error_halfedge, m_pRemesh) != NULL) {
      if (track_information) {
        std::cout << "1 edge flipped" << std::endl;
      }
      return;
    }
    // step 2: try to relocate
    int nb_relocate = relocate_applied(input_facet_tree, -1.0, max_error,
      -1.0, reduce_complexity, large_error_queue, small_value_queue,
      max_error_halfedge, m_pRemesh);
    if (nb_relocate > 0) {
      if (track_information) {
        std::cout << nb_relocate << " edges relocated" << std::endl;
      }
      return;
    }
    // step 3: try to split
    split_edge(input_facet_tree, -1.0, max_error, -1.0, reduce_complexity,
      large_error_queue, small_value_queue, max_error_halfedge, m_pRemesh);
    if (track_information) {
      std::cout << "1 edge splitted" << std::endl;
    }
  }

  // split
  Vertex_handle split_edge(const Facet_tree &input_facet_tree,
      FT max_error_threshold, FT max_error, FT min_radian,
      bool reduce_complexity, DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue,
      Halfedge_handle hh, Polyhedron *m_pRemesh) const {
    // max_error > 0 means we want to reduce error; othewiese improve radian
    // step 1: backup the original in_links and the edge types
    std::set<Facet_handle> one_ring_facets, extended_facets;
    one_ring_facets.insert(hh->facet());
    if (!hh->opposite()->is_border()) {
      one_ring_facets.insert(hh->opposite()->facet());
    }
    m_pRemesh->extend_facets(one_ring_facets, m_stencil_ring_size,
      &extended_facets);
    Link_iter_list facet_in_links, edge_in_links;
    Link_pointer_list vertex_in_links;
    backup_local_in_links(extended_facets, &facet_in_links, &edge_in_links,
      &vertex_in_links);
    std::map<Vertex_handle, bool> crease_map;
    if (m_inherit_element_types) {
      bool is_crease = hh->normal_dihedral() == -1.0 ?
        hh->opposite()->is_crease() : hh->is_crease();
      crease_map[m_pRemesh->get_source_vertex(hh)] = is_crease;
      crease_map[m_pRemesh->get_target_vertex(hh)] = is_crease;
    }
    // step 2: remove from queue if necessary
    if (small_value_queue != NULL) {
      remove_small_value_edges_before_split(hh, large_error_queue,
        small_value_queue, m_pRemesh);
    }
    // step 3: split the edge (also update edge types)
    Point new_point = m_pRemesh->midpoint(hh);
    Halfedge_handle hnew = m_pRemesh->split_long_edge(new_point, hh);
    Vertex_handle vh = hnew->vertex();
    if (m_inherit_element_types) {
      Halfedge_around_vertex_circulator vcirc = vh->vertex_begin();
      Halfedge_around_vertex_circulator vend = vcirc;
      typename std::map<Vertex_handle, bool>::iterator it;
      CGAL_For_all(vcirc, vend) {
        it = crease_map.find(m_pRemesh->get_source_vertex(vcirc));
        if (it != crease_map.end()) {   // splited halfedge
          Halfedge_handle hh = vcirc;
          if (hh->is_border()) {
            hh = hh->opposite();
          }
          hh->is_crease() = it->second;
          hh->normal_dihedral() = -0.5;               // temperary mark
          hh->opposite()->normal_dihedral() = -1.0;   // make opposite empty
        }
      }
    }
    // step 4: update local links (optimize the position if necessary)
    one_ring_facets.clear();
    extended_facets.clear();
    m_pRemesh->collect_one_ring_facets_incident_to_vertex(vh,
      &one_ring_facets);
    m_pRemesh->extend_facets(one_ring_facets, m_stencil_ring_size,
      &extended_facets);
    generate_local_links(input_facet_tree, false,
      &facet_in_links, &edge_in_links, &vertex_in_links,
      vh, &extended_facets, m_pRemesh);
    if (max_error > 0 || m_optimize_after_local_operations) {
      // if reduce error, we definitely optimize; otherweise it depends
      optimize_vertex_position(input_facet_tree, &facet_in_links,
        &edge_in_links, &vertex_in_links, vh,
        &extended_facets, m_pRemesh);
    }
    m_pRemesh->calculate_max_squared_errors(&extended_facets);
    // step 5: add to queue if necessary
    if (small_value_queue != NULL) {
      add_small_value_edges_after_split(vh, max_error_threshold,
        reduce_complexity, large_error_queue, small_value_queue, m_pRemesh);
    }
    // step 6: update normals
    std::set<Vertex_handle> one_ring_vertices;
    m_pRemesh->collect_vertices(one_ring_facets, &one_ring_vertices);
    m_pRemesh->calculate_local_normals(&one_ring_facets, &one_ring_vertices);
    // step 7: flip or relocate if necessary
    if (max_error < 0) {
      if (m_flip_after_split_and_collapse) {
        Halfedge_list halfedges;
        Halfedge_around_vertex_circulator vcirc = vh->vertex_begin();
        Halfedge_around_vertex_circulator vend = vcirc;
        CGAL_For_all(vcirc, vend) {
          if (!vcirc->is_border()) {
            halfedges.push_back(vcirc->prev());
          }
        }
        flip_edges(input_facet_tree, max_error_threshold, max_error,
          min_radian, reduce_complexity, large_error_queue,
          small_value_queue, &halfedges, m_pRemesh);
      }
      else if (m_relocate_after_local_operations) {
        std::set<Vertex_handle> one_ring_vertices;
        m_pRemesh->collect_indicent_vertices(vh, &one_ring_vertices);
        relocate_vertices(input_facet_tree, max_error_threshold, max_error,
          min_radian, reduce_complexity, large_error_queue,
          small_value_queue, &one_ring_vertices, m_pRemesh);
      }
    }
    return vh;
  }

  // collapse
  Vertex_handle collapse_applied(const Facet_tree &input_facet_tree,
      FT max_error_threshold, FT min_radian, bool reduce_complexity,
      bool *infinite_loop, DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue, Halfedge_handle hh,
      Polyhedron *m_pRemesh) {
    /* if min_radian > 0, we improve min_radian;
    if infinite_loop is not NULL, we check the infinite loop case;
    if improve_min_radian, we improve the min radian; otherwise,
    we collapse to reduce the mesh complexity */
    // step 1: topology constraints check
    if (!m_pRemesh->is_collapsible(hh)) {
      return NULL;
    }
    // step 2: geometry constraints check
    Halfedge_list halfedges;  // use the halfedges to represent facets
    bool is_ring = m_pRemesh->predict_facets_after_collapse(hh, &halfedges);
    Point new_point = get_initial_point_for_collapse(
      m_inherit_element_types, m_feature_control_delta, hh, m_pRemesh);
    if (m_keep_vertex_in_one_ring &&
      m_pRemesh->collapse_would_cause_wrinkle(halfedges, new_point, hh)) {
      return NULL;
    }
    // step 3: backup the original local links
    std::set<Facet_handle> one_ring_facets, extended_facets;
    m_pRemesh->collect_one_ring_facets_incident_to_edge(hh, &one_ring_facets);
    m_pRemesh->extend_facets(one_ring_facets, m_stencil_ring_size,
      &extended_facets);
    Link_iter_list facet_in_links, edge_in_links;
    Link_pointer_list vertex_in_links;
    Point_list facet_in_end_points, edge_in_end_points, vertex_in_end_points;
    backup_local_in_links(extended_facets, &facet_in_links,
      &facet_in_end_points, &edge_in_links, &edge_in_end_points,
      &vertex_in_links, &vertex_in_end_points);
    // step 4: simulate the edge collapse
    FT error = DOUBLE_MAX, radian = 0.0;
    simulate_edge_collapse(input_facet_tree, one_ring_facets,
      extended_facets, &halfedges, hh, is_ring,
      &facet_in_links, &edge_in_links, &vertex_in_links,
      &error, &radian, &new_point);
    restore_local_links(facet_in_end_points, &facet_in_links,
      edge_in_end_points, &edge_in_links,
      vertex_in_end_points, &vertex_in_links);
    // step 5: fidelity constraints check (max_error)
    if (error >= max_error_threshold) {
      return NULL;
    }
    // step 6: quality constraints check (min_radian) if necessary
    if (min_radian > 0 && radian < min_radian) {
      return NULL;
    }
    // step 7: infinite loops case check if necessary
    if (infinite_loop != NULL) {
      *infinite_loop = caused_infinite_loop(hh, m_pRemesh);
      if (*infinite_loop) {
        return NULL;
      }
    }
    // step 8: collapse the edge authentically
    Vertex_handle vh = collapse_edge(input_facet_tree, max_error_threshold,
      min_radian, reduce_complexity, large_error_queue, small_value_queue,
      &facet_in_links, &edge_in_links, &vertex_in_links, hh, new_point,
      m_pRemesh);
    return vh;
  }

  Point get_initial_point_for_collapse(bool inherit_element_types,
      FT feature_control_delta, Halfedge_handle hh,
      Polyhedron *m_pRemesh) const {
    // If inherit_element_type is true, we use the vertex_type;
    // Otherwise, we use the qem to get intial point for collapse
    if (inherit_element_types) {  // use the inherited vertex types
      Vertex_const_handle vp = m_pRemesh->get_source_vertex(hh);
      Vertex_const_handle vq = m_pRemesh->get_target_vertex(hh);
      Halfedge_const_list effective_edges;
      VertexType vtp = m_pRemesh->get_vertex_type(
        true, feature_control_delta, vp, &effective_edges);
      VertexType vtq = m_pRemesh->get_vertex_type(
        true, feature_control_delta, vq, &effective_edges);
      if (vtp == vtq) {
        // option: use dht qem instead
        return CGAL::midpoint(vp->point(), vq->point());
      }
      else {
        if (vtp == VertexType::k_feature_vertex) {
          return vp->point();
        }
        else if (vtp == VertexType::k_crease_vertex) {
          return vtq == VertexType::k_feature_vertex ? vq->point() :
            vp->point();
        }
        else {
          return vq->point();
        }
      }
    }
    else {                        // use the qem instead
      return m_pRemesh->get_least_qem_point(hh);
    }
    // backup option:
    /*FT fi_p = vp->feature_intensity();
    FT fi_q = vq->feature_intensity();
    if (CGAL::abs(fi_p - fi_q) <
    CGAL::max(fi_p, fi_q) * m_feature_difference_delta) {
    return midpoint(hh);
    }
    else {
    return fi_p > fi_q ? vp->point() : vq->point();
    }*/
  }

  void simulate_edge_collapse(const Facet_tree &input_facet_tree,
      const std::set<Facet_handle> &one_ring_facets,
      const std::set<Facet_handle> &extended_facets,
      Halfedge_list *halfedges, Halfedge_handle hh, bool is_ring,
      Link_iter_list *facet_in_links, Link_iter_list *edge_in_links,
      Link_pointer_list *vertex_in_links, FT *error, FT *radian,
      Point *new_point) const {
    // step 1: construct the local Remesh
    Polyhedron local_mesh;
    Vertex_handle vh = construct_local_mesh(one_ring_facets,
      extended_facets, halfedges, *new_point, is_ring, &local_mesh);
    // step 2: get the in link_facets (for function compatability)
    std::set<Facet_handle> in_link_facets;
    for (Facet_iterator fi = local_mesh.facets_begin();
      fi != local_mesh.facets_end(); ++fi) {
      in_link_facets.insert(fi);
    }
    generate_local_links(input_facet_tree, true,
      facet_in_links, edge_in_links, vertex_in_links,
      vh, &in_link_facets, &local_mesh);
    // step 3: optimize the vertex position if necessary
    if (m_optimize_after_local_operations) {
      optimize_vertex_position(input_facet_tree, facet_in_links,
        edge_in_links, vertex_in_links, vh, &in_link_facets, &local_mesh);
    }
    // step 4: update the max_errors for facets
    local_mesh.calculate_max_squared_errors(&in_link_facets);
    // step 5: calcualte the error, radian and new_point
    local_mesh.get_local_maximal_error(in_link_facets, error);
    *radian = local_mesh.get_minimal_radian_around_vertex(vh);
    *new_point = vh->point();
  }

  Vertex_handle construct_local_mesh(
      const std::set<Facet_handle> &one_ring_facets,
      const std::set<Facet_handle> &extended_facets,
      Halfedge_list *halfedges, const Point &new_point, bool is_ring,
      Polyhedron *local_mesh) const {
    // constructed the 2-manifold local mesh, and copy the feature intensities
    // step 1: construct the map (point -> Vertex_handle) and one_ring_facets
    std::map<Point, Vertex_handle> points_map;
    for (Halfedge_iter hi = halfedges->begin();
      hi != halfedges->end(); ++hi) {
      Halfedge_handle hh = *hi;
      points_map[hh->vertex()->point()] = hh->vertex();
    }
    if (!is_ring) {
      Halfedge_handle hh = *(halfedges->begin());
      points_map[hh->vertex()->point()] = hh->vertex();
    }
    // step 2: construct the manifold local mesh
    Halfedge_iter hi = halfedges->begin();
    Halfedge_handle hh = *hi;
    Halfedge_handle h1 = local_mesh->make_triangle(new_point,
      local_mesh->get_source_vertex(hh)->point(),
      local_mesh->get_target_vertex(hh)->point());
    Halfedge_handle h3 = h1->next()->opposite();
    h1 = h1->opposite();
    Halfedge_handle h2;
    ++hi;
    for (; hi != halfedges->end(); ++hi) {
      hh = *hi;
      h2 = local_mesh->make_triangle(new_point,
        local_mesh->get_source_vertex(hh)->point(),
        local_mesh->get_target_vertex(hh)->point());
      h2 = h2->next()->opposite();
      h1 = local_mesh->weld_edge(h1, h2);
      h1 = h1->opposite()->prev()->opposite();
    }
    if (is_ring) {
      h1 = local_mesh->weld_edge(h1, h3);
    }
    else {
      h1 = h1->opposite();
    }
    // step 3: copy the old feature intensities and calcualte the new one
    Halfedge_around_vertex_circulator vcirc = h1->vertex()->vertex_begin();
    Halfedge_around_vertex_circulator vend = vcirc;
    CGAL_For_all(vcirc, vend) {
      Vertex_handle vh = vcirc->opposite()->vertex();
      Vertex_handle old_vh = points_map[vh->point()];
      vh->gaussian_curvature() = old_vh->gaussian_curvature();
      vh->max_halfedge_dihedral() = old_vh->max_halfedge_dihedral();
    }
    // step 4: extend the local mesh
    // Now only the one_ring facets are manifold, others are triangle-soups
    std::set<Facet_handle> differ_facets;
    std::set_difference(extended_facets.begin(), extended_facets.end(),
      one_ring_facets.begin(), one_ring_facets.end(),
      std::inserter(differ_facets, differ_facets.end()));
    for (typename std::set<Facet_handle>::iterator it = differ_facets.begin();
      it != differ_facets.end(); ++it) {
      Facet_handle fh = *it;
      hh = fh->halfedge();
      local_mesh->make_triangle(hh->vertex()->point(),
        hh->next()->vertex()->point(),
        hh->prev()->vertex()->point());
    }
    return h1->vertex();
  }

  bool caused_infinite_loop(Halfedge_handle hh, Polyhedron *m_pRemesh) {
    typename std::map<Point, std::map<FT, Visit_iter>, Point_Comp>::iterator it1;
    Visit_iter it2;
    Point point = m_pRemesh->get_opposite_vertex(hh)->point();
    FT sl = m_pRemesh->squared_length(hh);
    sl = to_approximation(sl);
    bool found = false;
    it1 = m_collapsed_map.find(point);
    if (it1 != m_collapsed_map.end()) {
      typename std::map<FT, Visit_iter>::iterator it = it1->second.find(sl);
      if (it != it1->second.end()) {
        found = true;
        it2 = it->second;
      }
    }
    if (found) {
      if (m_track_information) {
        std::cout << "Point(" << point << ") with length "
          << sl << ": collapse denied.";
      }
      //option: we may delete it instead of splicing front
      m_collapsed_list.splice(m_collapsed_list.begin(),
        m_collapsed_list, it2);
      return true;
    }
    else {
      m_collapsed_list.push_front(std::pair<Point, FT>(point, sl));
      m_collapsed_map[point][sl] = m_collapsed_list.begin();
      while (m_collapsed_list.size() > m_collapsed_list_size) {
        const Point &p = m_collapsed_list.back().first;
        FT sl = m_collapsed_list.back().second;
        it1 = m_collapsed_map.find(p);
        if (it1 != m_collapsed_map.end()) {
          it1->second.erase(sl);
          if (it1->second.size() == 0) {
            m_collapsed_map.erase(it1);
          }
        }
        m_collapsed_list.pop_back();
      }
      return false;
    }
  }

  Vertex_handle collapse_edge(const Facet_tree &input_facet_tree,
      FT max_error_threshold, FT min_radian, bool reduce_complexity,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue,
      Link_iter_list *facet_in_links, Link_iter_list *edge_in_links,
      Link_pointer_list *vertex_in_links, Halfedge_handle hh,
      Point new_point, Polyhedron *m_pRemesh) const {
    // step 1: backup the edge type if necessary
    std::map<Vertex_handle, bool> crease_map;
    if (m_inherit_element_types) {
      Vertex_handle vh = m_pRemesh->get_opposite_vertex(hh);
      crease_map[vh] = m_pRemesh->is_crease_edge(hh->prev()) ||
        m_pRemesh->is_crease_edge(hh->next());
      if (!hh->opposite()->is_border()) {
        vh = m_pRemesh->get_opposite_vertex(hh->opposite());
        crease_map[vh] = m_pRemesh->is_crease_edge(hh->opposite()->prev()) ||
          m_pRemesh->is_crease_edge(hh->opposite()->next());
      }
    }
    // step 2: remove from queue if necessary
    if (small_value_queue != NULL) {
      remove_small_value_edges_before_collapse(hh, large_error_queue,
        small_value_queue, m_pRemesh);
    }
    // step 3: collapse the edge (also update edge types here)
    Vertex_handle v_joined = m_pRemesh->collapse_short_edge(new_point, hh);
    if (v_joined == NULL) {
      return NULL;
    }
    if (m_inherit_element_types) {
      Halfedge_around_vertex_circulator vcirc = v_joined->vertex_begin();
      Halfedge_around_vertex_circulator vend = vcirc;
      typename std::map<Vertex_handle, bool>::iterator it;
      CGAL_For_all(vcirc, vend) {
        it = crease_map.find(m_pRemesh->get_source_vertex(vcirc));
        if (it != crease_map.end()) { // the merged halfedge
          Halfedge_handle hh = vcirc;
          if (hh->is_border()) {
            hh = hh->opposite();
          }
          hh->is_crease() = it->second;
          hh->normal_dihedral() = -0.5;             // temperary mark
          hh->opposite()->normal_dihedral() = -1.0; // mark opposite as empty
        }
      }
    }
    // step 4: update local links (optimize already performed in simulation)
    std::set<Facet_handle> one_ring_facets, extended_facets;
    m_pRemesh->collect_one_ring_facets_incident_to_vertex(v_joined,
      &one_ring_facets);
    m_pRemesh->extend_facets(one_ring_facets, m_stencil_ring_size,
      &extended_facets);
    generate_local_links(input_facet_tree, true,
      facet_in_links, edge_in_links, vertex_in_links,
      v_joined, &extended_facets, m_pRemesh);
    m_pRemesh->calculate_max_squared_errors(&extended_facets);
    // step 5: add to queue if necessary
    if (small_value_queue != NULL) {
      add_small_value_edges_after_collapse(v_joined, max_error_threshold,
        reduce_complexity, large_error_queue, small_value_queue, m_pRemesh);
    }
    // step 6: update normals
    std::set<Vertex_handle> one_ring_vertices;
    m_pRemesh->collect_vertices(one_ring_facets, &one_ring_vertices);
    m_pRemesh->calculate_local_normals(&one_ring_facets, &one_ring_vertices);

    // step 7: flip or relocate if necessary
    if (min_radian > 0) {
      if (m_flip_after_split_and_collapse) {
        Halfedge_list halfedges;
        Halfedge_around_vertex_circulator vcirc = v_joined->vertex_begin();
        Halfedge_around_vertex_circulator vend = vcirc;
        CGAL_For_all(vcirc, vend) {
          if (!vcirc->is_border()) {
            halfedges.push_back(vcirc);
          }
        }
        flip_edges(input_facet_tree, max_error_threshold, -1.0,
          min_radian, reduce_complexity, large_error_queue,
          small_value_queue, &halfedges, m_pRemesh);
      }
      else if (m_relocate_after_local_operations) {
        std::set<Vertex_handle> one_ring_vertices;
        m_pRemesh->collect_indicent_vertices(v_joined, &one_ring_vertices);
        relocate_vertices(input_facet_tree, max_error_threshold, -1.0,
          min_radian, reduce_complexity, large_error_queue,
          small_value_queue, &one_ring_vertices, m_pRemesh);
      }
    }
    return v_joined;
  }

  // flip
  int flip_applied(const Facet_tree &input_facet_tree,
      FT max_error_threshold, FT max_error, FT min_radian,
      bool reduce_complexity, DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_radian_queue, Halfedge_handle hh,
      Polyhedron *m_pRemesh) const {
    // step 1: construct the incident_edges
    Halfedge_list incident_edges;   // two incident edges will be included
    FT sd_prev = m_pRemesh->squared_length(hh->prev());
    FT sd_next = m_pRemesh->squared_length(hh->next());
    if (sd_prev >= sd_next) {
      incident_edges.push_back(hh->prev());
      incident_edges.push_back(hh->next());
    }
    else {
      incident_edges.push_back(hh->next());
      incident_edges.push_back(hh->prev());
    }
    // step 2: try to flip
    return flip_edges(input_facet_tree, max_error_threshold, max_error,
      min_radian, reduce_complexity, large_error_queue, small_radian_queue,
      &incident_edges, m_pRemesh);
  }

  int flip_edges(const Facet_tree &input_facet_tree, FT max_error_threshold,
      FT max_error, FT min_radian, bool reduce_complexity,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_radian_queue, Halfedge_list *halfedges,
      Polyhedron *m_pRemesh) const {
    int num = 0;
    for (Halfedge_iter it = halfedges->begin();
      it != halfedges->end(); ++it) {
      Halfedge_handle hh = *it;
      Halfedge_handle hnew = flip_edge(input_facet_tree, max_error_threshold,
        max_error, min_radian, reduce_complexity, large_error_queue,
        small_radian_queue, hh, m_pRemesh);
      if (hnew != NULL) {
        ++num;
      }
    }
    return num;
  }

  Halfedge_handle flip_edge(const Facet_tree &input_facet_tree,
      FT max_error_threshold, FT max_error, FT min_radian,
      bool reduce_complexity, DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue, Halfedge_handle hh,
      Polyhedron *m_pRemesh) const {
    // max_error > 0 means we want to reduce error; otherwise improve radian
    // step 1: topology constraints and min_radian constraints check
    if (!m_pRemesh->is_flippable(m_inherit_element_types,
      m_edge_flip_strategy, hh)) {
      return NULL;
    }
    // step 2: backup the original local links
    std::set<Facet_handle> one_ring_facets, extended_facets;
    one_ring_facets.insert(hh->facet());
    one_ring_facets.insert(hh->opposite()->facet());
    m_pRemesh->extend_facets(one_ring_facets, m_stencil_ring_size,
      &extended_facets);
    Link_iter_list facet_in_links, edge_in_links;
    Link_pointer_list vertex_in_links;
    backup_local_in_links(extended_facets, &facet_in_links,
      &edge_in_links, &vertex_in_links);
    // step 3: remove from queue if necessary
    if (small_value_queue != NULL) {
      remove_small_value_edges_before_flip(hh, large_error_queue,
        small_value_queue, m_pRemesh);
    }
    // step 4: flip (do not optimize the vertex position)
    Halfedge_handle hnew = m_pRemesh->flip_edge(hh);
    // step 5: update local links
    one_ring_facets.clear();
    extended_facets.clear();
    one_ring_facets.insert(hnew->facet());
    one_ring_facets.insert(hnew->opposite()->facet());
    m_pRemesh->extend_facets(one_ring_facets, m_stencil_ring_size,
      &extended_facets);
    generate_local_links(input_facet_tree, false, &facet_in_links,
      &edge_in_links, &vertex_in_links, hnew, &extended_facets, m_pRemesh);
    m_pRemesh->calculate_max_squared_errors(&extended_facets);
    FT error_after;
    m_pRemesh->get_local_maximal_error(extended_facets, &error_after);
    // step 6: rollback if the flip violated the max_error constraints;
    bool flipped;
    if (max_error > 0) {
      flipped = error_after < max_error;
    }
    else {
      flipped = error_after < max_error_threshold;
    }
    if (!flipped) {
      hh = m_pRemesh->flip_edge(hnew);
      one_ring_facets.clear();
      extended_facets.clear();
      one_ring_facets.insert(hh->facet());
      one_ring_facets.insert(hh->opposite()->facet());
      m_pRemesh->extend_facets(one_ring_facets, m_stencil_ring_size,
        &extended_facets);
      generate_local_links(input_facet_tree, false, &facet_in_links,
        &edge_in_links, &vertex_in_links, hh, &extended_facets, m_pRemesh);
      m_pRemesh->calculate_max_squared_errors(&extended_facets);
    }
    // step 7: add to queue if necessary
    if (small_value_queue != NULL) {
      if (flipped) {
        add_small_value_edges_after_flip(hnew, max_error_threshold,
          reduce_complexity, large_error_queue, small_value_queue,
          m_pRemesh);
      }
      else {
        add_small_value_edges_after_flip(hh, max_error_threshold,
          reduce_complexity, large_error_queue, small_value_queue,
          m_pRemesh);
      }
    }
    // step 8: update the normals
    std::set<Vertex_handle> one_ring_vertices;
    m_pRemesh->collect_vertices(one_ring_facets, &one_ring_vertices);
    m_pRemesh->calculate_local_normals(&one_ring_facets, &one_ring_vertices);
    // step 9: relocate if necessary (only when we want to improve angle)
    if (max_error < 0) {
      if (flipped && m_relocate_after_local_operations) {
        std::set<Vertex_handle> vertices;
        vertices.insert(m_pRemesh->get_source_vertex(hnew));
        vertices.insert(m_pRemesh->get_target_vertex(hnew));
        vertices.insert(m_pRemesh->get_opposite_vertex(hnew));
        vertices.insert(m_pRemesh->get_opposite_vertex(hnew->opposite()));
        std::set<Facet_handle> facets;  // used to calculate local_min_radian
        for (auto it = vertices.begin(); it != vertices.end(); ++it) {
          m_pRemesh->collect_one_ring_facets_incident_to_vertex(*it, &facets);
        }
        FT local_min_radian = m_pRemesh->get_local_minimal_radian(facets);
        relocate_vertices(input_facet_tree, max_error_threshold, max_error,
          local_min_radian, reduce_complexity, large_error_queue,
          small_value_queue, &vertices, m_pRemesh);
      }
    }
    return flipped ? hnew : NULL;
  }

  // relocate
  int relocate_applied(const Facet_tree &input_facet_tree,
      FT max_error_threshold, FT max_error, FT min_radian,
      bool reduce_complexity, DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue,
      Halfedge_handle hh, Polyhedron *m_pRemesh) const {
    // If max_error > 0, we reduce error; otherwise, we improve min_radian
    // If reduce_complexity, we update collapse_candidate_queue;
    //    otherwise, we update small_radian_queue.
    // step 1: construct the vertices
    std::set<Vertex_handle> incident_vertices;
    incident_vertices.insert(m_pRemesh->get_opposite_vertex(hh));
    incident_vertices.insert(m_pRemesh->get_target_vertex(hh));
    incident_vertices.insert(m_pRemesh->get_source_vertex(hh));
    // step 2: relocate these incident vertices
    return relocate_vertices(input_facet_tree, max_error_threshold,
      max_error, min_radian, reduce_complexity, large_error_queue,
      small_value_queue, &incident_vertices, m_pRemesh);
  }

  int relocate_vertices(const Facet_tree &input_facet_tree,
      FT max_error_threshold, FT max_error, FT min_radian,
      bool reduce_complexity, DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_radian_queue,
      std::set<Vertex_handle> *vertices, Polyhedron *m_pRemesh) const {
    // step 1: construct the vertex_map to make bigger distance relocate first
    std::map<FT, Vertex_handle> vertex_map;
    for (auto it = vertices->begin(); it != vertices->end(); ++it) {
      Vertex_handle vh = *it;
      Point p = get_initial_point_for_relocate(input_facet_tree,
        m_inherit_element_types, m_max_samples_per_area,
        m_feature_control_delta, vh, m_pRemesh);
      vertex_map[CGAL::squared_distance(vh->point(), p)] = vh;
    }
    // step 2: relocate these vertices in order
    int num = 0;
    for (auto it = vertex_map.rbegin(); it != vertex_map.rend(); ++it) {
      Vertex_handle vh = it->second;
      Point initial_point = get_initial_point_for_relocate(input_facet_tree,
        m_inherit_element_types, m_max_samples_per_area,
        m_feature_control_delta, vh, m_pRemesh);
      num += relocate_vertex(input_facet_tree, max_error_threshold, max_error,
        min_radian, reduce_complexity, large_error_queue, small_radian_queue,
        initial_point, vh, m_pRemesh);
    }
    return num;
  }

  Point get_initial_point_for_relocate(const Facet_tree &facet_tree,
      bool inherit_element_types, FT max_samples_per_area,
      FT feature_control_delta, Vertex_handle vh,
      Polyhedron *m_pRemesh) const {
    FT sum_area = 0.0;
    Halfedge_around_vertex_const_circulator vcirc = vh->vertex_begin();
    Halfedge_around_vertex_const_circulator vend = vcirc;
    CGAL_For_all(vcirc, vend) {
      if (!vcirc->is_border()) {
        sum_area += m_pRemesh->area(vcirc->facet());
      }
    }
    Point initial_point;
    if (sum_area <= 1.0 / max_samples_per_area * vh->degree()) {
      // if the area is too small, return its average center
      initial_point = m_pRemesh->average_center(vh);
    }
    else {
      if (m_relocate_strategy == RelocateStrategy::k_cvt_barycenter) {
        initial_point = m_pRemesh->cvt_barycenter(inherit_element_types,
          feature_control_delta, vh);
      }
      else {
        initial_point = m_pRemesh->barycenter(inherit_element_types,
          feature_control_delta, vh);
      }
    }
    // if it is a on a crease, then it is already on the facet_tree
    return facet_tree.closest_point(initial_point);
  }

  bool relocate_vertex(const Facet_tree &input_facet_tree,
      FT max_error_threshold, FT max_error, FT min_radian,
      bool reduce_complexity, DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue, const Point &new_point,
      Vertex_handle vh, Polyhedron *m_pRemesh) const {
    // If max_error > 0, we want to reduce error; otherwiese increase radian
    // If reduce_complexity, we update collapse_candidate_queue;
    //    otherwise, we update small_radian_queue.
    // step 1: geometry constraints check
    if (CGAL::squared_distance(vh->point(), new_point) < SQUARED_MIN_VALUE) {
      return false;
    }
    if (m_keep_vertex_in_one_ring &&
      m_pRemesh->relocate_would_cause_wrinkle(new_point, vh)) {
      return false;
    }
    // step 2: backup the original local links
    std::set<Facet_handle> one_ring_facets, extended_facets;
    m_pRemesh->collect_one_ring_facets_incident_to_vertex(vh,
      &one_ring_facets);
    m_pRemesh->extend_facets(one_ring_facets, m_stencil_ring_size,
      &extended_facets);
    Link_iter_list facet_in_links, edge_in_links;
    Link_pointer_list vertex_in_links;
    std::map<Facet_handle, Link_list> facet_out_map;
    std::map<Halfedge_handle, Link_list> edge_out_map;
    Link vertex_out;
    backup_local_links(extended_facets, &facet_out_map, &edge_out_map,
      &vertex_out, &facet_in_links, &edge_in_links, &vertex_in_links, vh);
    // step 3: remove from queue if necessary
    if (small_value_queue != NULL) {
      remove_small_value_edges_before_relocate(vh, large_error_queue,
        small_value_queue, m_pRemesh);
    }
    // step 4: relocate the vertex position
    Point old_point = vh->point(); // backup position in case rolling back
    vh->point() = new_point;
    // step 5: update local links (optimized included if necessary)
    generate_local_links(input_facet_tree, false,
      &facet_in_links, &edge_in_links, &vertex_in_links,
      vh, &extended_facets, m_pRemesh);
    if (max_error > 0 || m_optimize_after_local_operations) {
      // if reduce error, we definitely optimize; otherweise it depends
      optimize_vertex_position(input_facet_tree, &facet_in_links,
        &edge_in_links, &vertex_in_links, vh, &extended_facets, m_pRemesh);
    }
    m_pRemesh->calculate_max_squared_errors(&extended_facets);
    FT error, radian;
    m_pRemesh->get_local_maximal_error(extended_facets, &error);
    radian = m_pRemesh->get_minimal_radian_around_vertex(vh);
    // step 6: rollback if the max_error or min_radian is violated
    FT radian_delta = m_pRemesh->to_radian(m_smooth_angle_delta);
    bool relocated;
    if (max_error > 0) {      // only reduce the error
      relocated = (error < max_error);
    }
    else {                    // improve the min_radian
      relocated = (error < max_error_threshold &&
        radian >= min_radian + radian_delta);
    }
    if (!relocated) {
      vh->point() = old_point;
      restore_local_links(facet_out_map, edge_out_map, vertex_out,
        &facet_in_links, &edge_in_links, &vertex_in_links,
        vh, &extended_facets, m_pRemesh);
      m_pRemesh->calculate_max_squared_errors(&extended_facets);
    }
    // step 7: add the queue if necesary
    if (small_value_queue != NULL) {
      add_small_value_edges_after_relocate(vh, max_error_threshold,
        reduce_complexity, large_error_queue, small_value_queue, m_pRemesh);
    }
    // step 8: update the normals
    std::set<Vertex_handle> one_ring_vertices;
    m_pRemesh->collect_vertices(one_ring_facets, &one_ring_vertices);
    m_pRemesh->calculate_local_normals(&one_ring_facets, &one_ring_vertices);

    return relocated;
  }

  // optimize
  void optimize_vertex_position(const Facet_tree &input_facet_tree,
      Link_iter_list *facet_in_links, Link_iter_list *edge_in_links,
      Link_pointer_list *vertex_in_links, Vertex_handle vh,
      std::set<Facet_handle> *in_link_facets, Polyhedron *m_pRemesh) const {
    // precondition: the samples are generated
    for (int i = 0; i < m_vertex_optimize_count; ++i) {
      // step 1: optimize the position of vh
      if (!optimize_vertex_by_one_ring(input_facet_tree, vh, m_pRemesh)) {
        break;
      }
      // step 2: update the samples
      generate_local_links(input_facet_tree, false, facet_in_links,
        edge_in_links, vertex_in_links, vh, in_link_facets, m_pRemesh);
    }
  }

  bool optimize_vertex_by_one_ring(const Facet_tree &input_facet_tree,
      Vertex_handle vh, Polyhedron *m_pRemesh) const {
    // step 1: check whether it is too small to optimize
    size_t nb_out_links = 0, nb_in_links = 0;
    Halfedge_around_vertex_circulator vcirc = vh->vertex_begin();
    Halfedge_around_vertex_circulator vend = vcirc;
    CGAL_For_all(vcirc, vend) {
      if (!vcirc->is_border()) {
        Facet_handle fh = vcirc->facet();
        nb_out_links += fh->facet_out_links().size();
        nb_in_links += fh->facet_in_links().size();
        nb_in_links += fh->edge_in_links().size();
        nb_in_links += fh->vertex_in_links().size();
      }
    }
    if (nb_out_links == 0 || nb_in_links == 0) {
      vh->point() = m_pRemesh->barycenter(m_inherit_element_types,
        m_feature_control_delta, vh);
      return false;     // the area is too small
    }
    // step 2: calcualte the optimized position
    Point new_point = get_optimized_position(vh, m_pRemesh);
    if (CGAL::squared_distance(new_point, vh->point()) < SQUARED_MIN_VALUE) {
      return false;     // the optimization is too little
    }
    if (m_keep_vertex_in_one_ring) {
      FT min_sd = m_pRemesh->ge_min_squared_distance_in_one_ring_facets(vh);
      FT sd = CGAL::squared_distance(new_point, vh->point());
      if (sd > min_sd) {
        vh->point() = vh->point() +
          (new_point - vh->point()) * CGAL::sqrt(min_sd) / CGAL::sqrt(sd);
      }
    }
    else {
      vh->point() = new_point;
    }
    if (m_optimize_strategy == OptimizeStrategy::k_Interpolation) {
      vh->point() = input_facet_tree.closest_point(vh->point());
    }
    return true;
  }

  Point get_optimized_position(Vertex_handle vh, Polyhedron *m_pRemesh) const {
    FT denominator = 0.0;
    Vector vec_sum = CGAL::NULL_VECTOR;
    bool use_facet_in_links = false, use_facet_out_links = false;
    bool use_edge_in_links = false, use_edge_out_links = false;
    bool use_vertex_in_links = false, use_vertex_out_links = false;
    // step 1: collect the facet samples
    use_facet_in_links =
      m_facet_optimize_type == OptimizeType::k_input_to_remesh ||
      m_facet_optimize_type == OptimizeType::k_both;
    use_facet_out_links =
      m_facet_optimize_type == OptimizeType::k_remesh_to_input ||
      m_facet_optimize_type == OptimizeType::k_both;
    // step 2: collect the edge samples
    use_edge_in_links =
      m_edge_optimize_type == OptimizeType::k_input_to_remesh ||
      m_edge_optimize_type == OptimizeType::k_both;
    use_edge_out_links =
      m_edge_optimize_type == OptimizeType::k_remesh_to_input ||
      m_edge_optimize_type == OptimizeType::k_both;
    // step 3: collect the vertex samples
    use_vertex_in_links =
      m_vertex_optimize_type == OptimizeType::k_input_to_remesh ||
      m_vertex_optimize_type == OptimizeType::k_both;
    use_vertex_out_links =
      m_vertex_optimize_type == OptimizeType::k_remesh_to_input ||
      m_vertex_optimize_type == OptimizeType::k_both;
    // step 4: get the optimized vector
    vec_sum = get_optimize_vector(use_facet_in_links, use_facet_out_links,
      use_edge_in_links, use_edge_out_links, use_vertex_in_links,
      use_vertex_out_links, vh, &denominator, m_pRemesh);
    // step 5: get the vertex position if valid
    if (denominator > SQUARED_MIN_VALUE) {
      vec_sum = vec_sum / denominator;
      Point new_point = CGAL::ORIGIN + vec_sum;
      return vh->point() +
        m_vertex_optimize_ratio * (new_point - vh->point());
    }
    else {
      return vh->point();
    }
  }

  Vector get_optimize_vector(bool use_facet_in_links,
      bool use_facet_out_links, bool use_edge_in_links,
      bool use_edge_out_links, bool use_vertex_in_links,
      bool use_vertex_out_links, Vertex_handle vh,
      FT *denominator, Polyhedron *m_pRemesh) const {
    Vector sum_vec = CGAL::NULL_VECTOR;
    FT weight = 0.0;
    Halfedge_around_vertex_const_circulator vcirc = vh->vertex_begin();
    Halfedge_around_vertex_const_circulator vend = vcirc;
    CGAL_For_all(vcirc, vend) {
      // vertex out link
      if (use_vertex_out_links) {
        const Point &out_start = vh->vertex_out_link().second.first;
        const Point &out_end = vh->vertex_out_link().second.second;
        weight = CGAL::sqrt(CGAL::squared_distance(out_start, out_end));
        if (m_use_feature_intensity_weights) {
          weight *= vh->vertex_out_link().first;
        }
        sum_vec = sum_vec + weight * (out_end - CGAL::ORIGIN);
        *denominator += weight;
      }
      // edge out links
      if (use_edge_out_links) {
        const Point &v0 = m_pRemesh->get_target_vertex(vcirc)->point();
        const Point &v1 = m_pRemesh->get_source_vertex(vcirc)->point();
        Halfedge_const_handle hh = vcirc;
        if (hh->normal_dihedral() == -1.0) {
          hh = hh->opposite();
        }
        FT length = m_pRemesh->length(hh);
        if (length >= MIN_VALUE) {            // we ignore too short edges
          for (Link_list_const_iter it = hh->edge_out_links().begin();
            it != hh->edge_out_links().end(); ++it) {
            const Point &out_start = it->second.first;
            const Point &out_end = it->second.second;
            weight = CGAL::sqrt(CGAL::squared_distance(out_start, out_end));
            if (m_use_feature_intensity_weights) {
              weight *= it->first;
            }
            FT alpha =
              CGAL::sqrt(CGAL::squared_distance(out_start, v1)) / length;
            sum_vec = sum_vec + weight * alpha * (out_end - CGAL::ORIGIN)
              - weight * alpha * (1.0 - alpha) * (v1 - CGAL::ORIGIN);
            *denominator += weight * alpha * alpha;
          }
        }
      }
      if (vcirc->is_border()) {
        continue;
      }
      Facet_const_handle fh = vcirc->facet();
      const Point &v0 = m_pRemesh->get_target_vertex(vcirc)->point();
      const Point &v1 = m_pRemesh->get_opposite_vertex(vcirc)->point();
      const Point &v2 = m_pRemesh->get_source_vertex(vcirc)->point();
      const FT facet_area = m_pRemesh->area(fh);
      if (facet_area < SQUARED_MIN_VALUE) {   // we ignore too small facets
        continue;
      }
      // facet out links
      if (use_facet_out_links) {
        for (Link_list_const_iter it = fh->facet_out_links().begin();
          it != fh->facet_out_links().end(); ++it) {
          const Point &out_start = it->second.first;
          const Point &out_end = it->second.second;
          weight = CGAL::sqrt(CGAL::squared_distance(out_start, out_end));
          if (m_use_feature_intensity_weights) {
            weight *= it->first;
          }
          FT alpha_0 = m_pRemesh->area(v1, v2, out_start) / facet_area;
          FT alpha_1 = m_pRemesh->area(v2, v0, out_start) / facet_area;
          FT alpha_2 = m_pRemesh->area(v0, v1, out_start) / facet_area;
          sum_vec = sum_vec + weight * alpha_0 * (out_end - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_1 * (v1 - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_2 * (v2 - CGAL::ORIGIN);
          *denominator += weight * alpha_0 * alpha_0;
        }
      }
      // vertex in links
      if (use_vertex_in_links) {
        for (Link_pointer_const_iter it = fh->vertex_in_links().begin();
          it != fh->vertex_in_links().end(); ++it) {
          Link *link = *it;
          const Point &in_start = link->second.first;
          const Point &in_end = link->second.second;
          weight = CGAL::sqrt(CGAL::squared_distance(in_start, in_end));
          if (m_use_feature_intensity_weights) {
            weight *= link->first;
          }
          FT alpha_0 = m_pRemesh->area(v1, v2, in_end) / facet_area;
          FT alpha_1 = m_pRemesh->area(v2, v0, in_end) / facet_area;
          FT alpha_2 = m_pRemesh->area(v0, v1, in_end) / facet_area;
          sum_vec = sum_vec + weight * alpha_0 * (in_start - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_1 * (v1 - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_2 * (v2 - CGAL::ORIGIN);
          *denominator += weight * alpha_0 * alpha_0;
        }
      }
      // edge in links
      if (use_edge_in_links) {
        for (Link_iter_list_const_iter it = fh->edge_in_links().begin();
          it != fh->edge_in_links().end(); ++it) {
          Link_list_const_iter lit = *it;
          const Point &in_start = lit->second.first;
          const Point &in_end = lit->second.second;
          weight = CGAL::sqrt(CGAL::squared_distance(in_start, in_end));
          if (m_use_feature_intensity_weights) {
            weight *= lit->first;
          }
          FT alpha_0 = m_pRemesh->area(v1, v2, in_end) / facet_area;
          FT alpha_1 = m_pRemesh->area(v2, v0, in_end) / facet_area;
          FT alpha_2 = m_pRemesh->area(v0, v1, in_end) / facet_area;
          sum_vec = sum_vec + weight * alpha_0 * (in_start - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_1 * (v1 - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_2 * (v2 - CGAL::ORIGIN);
          *denominator += weight * alpha_0 * alpha_0;
        }
      }
      // facet in links
      if (use_facet_in_links) {
        for (Link_iter_list_const_iter it = fh->facet_in_links().begin();
          it != fh->facet_in_links().end(); ++it) {
          Link_list_const_iter lit = *it;
          const Point &in_start = lit->second.first;
          const Point &in_end = lit->second.second;
          weight = CGAL::sqrt(CGAL::squared_distance(in_start, in_end));
          if (m_use_feature_intensity_weights) {
            weight *= lit->first;
          }
          FT alpha_0 = m_pRemesh->area(v1, v2, in_end) / facet_area;
          FT alpha_1 = m_pRemesh->area(v2, v0, in_end) / facet_area;
          FT alpha_2 = m_pRemesh->area(v0, v1, in_end) / facet_area;
          sum_vec = sum_vec + weight * alpha_0 * (in_start - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_1 * (v1 - CGAL::ORIGIN)
            - weight * alpha_0 * alpha_2 * (v2 - CGAL::ORIGIN);
          *denominator += weight * alpha_0 * alpha_0;
        }
      }
    }
    return sum_vec;
  }

  // utilities
  inline FT to_approximation(FT value) const {
    FT precison = MAX_VALUE;
    int temp_value = value * precison;
    return temp_value / precison;
  }

  Facet_tree::Point_and_primitive_id get_closest_point_and_primitive(
      const std::set<Facet_handle> &in_link_facets,
      const Point &point, Polyhedron *polyhedron) const {
    // find the closest point in the facet set
    Facet_tree::Point_and_primitive_id pp;
    FT min_sd = DOUBLE_MAX;
    Point nearest_point;
    for (auto it = in_link_facets.begin(); it != in_link_facets.end(); ++it) {
      Facet_handle fh = *it;
      FT sd = polyhedron->squared_distance(point, fh, &nearest_point);
      if (sd < min_sd) {
        pp.first = nearest_point;
        pp.second = fh;
        min_sd = sd;
      }
    }
    return pp;
  }

 private:
  // general parameters
  FT m_max_error_threshold;
  FT m_min_angle_threshold;
  int m_max_mesh_complexity;
  FT m_smooth_angle_delta;
  bool m_apply_edge_flip;
  EdgeFlipStrategy m_edge_flip_strategy;
  bool m_flip_after_split_and_collapse;
  bool m_relocate_after_local_operations;
  RelocateStrategy m_relocate_strategy;
  bool m_keep_vertex_in_one_ring;
  bool m_use_local_aabb_tree;
  int m_collapsed_list_size;
  bool m_decrease_max_errors;
  bool m_track_information;
  bool m_apply_initial_mesh_simplification;
  bool m_apply_final_vertex_relocation;
  // sample parameters
  int m_samples_per_facet_in;
  int m_samples_per_facet_out;
  int m_max_samples_per_area;
  int m_min_samples_per_triangle;
  int m_bvd_iteration_count;
  SampleNumberStrategy m_sample_number_strategy;
  SampleStrategy m_sample_strategy;
  bool m_use_stratified_sampling;
  // feature function parameters
  FT m_sum_theta;
  FT m_sum_delta;
  FT m_dihedral_theta;
  FT m_dihedral_delta;
  FT m_feature_difference_delta;
  FT m_feature_control_delta;
  bool m_inherit_element_types;
  bool m_use_feature_intensity_weights;
  // vertex relocate_parameters
  int m_vertex_optimize_count;
  FT m_vertex_optimize_ratio;
  int m_stencil_ring_size;
  OptimizeStrategy m_optimize_strategy;
  OptimizeType m_facet_optimize_type;
  OptimizeType m_edge_optimize_type;
  OptimizeType m_vertex_optimize_type;
  bool m_optimize_after_local_operations;

  // the collapse operator
  Visit_list m_collapsed_list;
  std::map<Point, std::map<FT, Visit_iter>, Point_Comp> m_collapsed_map;
};

#endif // REMESH_H_

/*void optimize_halfedge_position(const Facet_tree &input_facet_tree,
  Link_iter_list *facet_in_links, Link_iter_list *edge_in_links,
  std::list<Link*> *vertex_in_links, Halfedge_handle hh,
  std::set<Facet_handle> *in_link_facets, Polyhedron *m_pRemesh) const {
  // precondition: the samples are generated
  for (int i = 0; i < m_vertex_optimize_count; ++i) {
    bool b1 = optimize_vertex_by_one_ring(input_facet_tree,
      local_mesh->get_source_vertex(hh), m_pRemesh);
    bool b2 = optimize_vertex_by_one_ring(input_facet_tree,
      local_mesh->get_target_vertex(hh), m_pRemesh);
    if (!b1 && !b2) {   // no vertex position optimized
      break;
    }
    generate_local_links(input_facet_tree, false, facet_in_links,
      edge_in_links, vertex_in_links, hh, in_link_facets, m_pRemesh);
  }
}
*/