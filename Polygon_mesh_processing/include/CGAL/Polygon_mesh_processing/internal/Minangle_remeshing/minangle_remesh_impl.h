// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Kaimo Hu

#ifndef SRC_INTERNAL_MINANGLE_REMESHING_MINANGLE_REMESH_IMPL_H_
#define SRC_INTERNAL_MINANGLE_REMESHING_MINANGLE_REMESH_IMPL_H_

// C/C++
#include <list>
#include <set>
#include <map>
#include <string>
#include <algorithm>
#include <utility>
// local
#include "mesh_properties.h"

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template<typename Kernel>
class Minangle_remesher {
 public:
  // type definitions
  // types
  typedef typename Mesh_properties<Kernel> Mesh_properties;
  typedef typename Mesh_properties::FT FT;
  typedef typename Mesh_properties::Vector Vector;
  typedef typename Mesh_properties::Normal Normal;
  typedef typename Mesh_properties::Point Point;
  typedef typename Mesh_properties::Point_Comp Point_Comp;
  typedef typename Mesh_properties::Point_pair Point_pair;
  typedef typename Mesh_properties::Bbox Bbox;
  typedef typename Mesh_properties::Mesh Mesh;
  // descriptors
  typedef typename Mesh_properties::halfedge_descriptor halfedge_descriptor;
  typedef typename Mesh_properties::edge_descriptor edge_descriptor;
  typedef typename Mesh_properties::vertex_descriptor vertex_descriptor;
  typedef typename Mesh_properties::face_descriptor face_descriptor;
  // element list
  typedef typename Mesh_properties::Halfedge_list Halfedge_list;
  typedef typename Mesh_properties::Halfedge_iter Halfedge_iter;
  typedef typename Mesh_properties::Halfedge_const_iter Halfedge_const_iter;
  typedef typename Mesh_properties::Edge_list Edge_list;
  typedef typename Mesh_properties::Edge_iter Edge_iter;
  typedef typename Mesh_properties::Edge_const_iter Edge_const_iter;
  typedef typename Mesh_properties::Vertex_list Vertex_list;
  typedef typename Mesh_properties::Vertex_iter Vertex_iter;
  typedef typename Mesh_properties::Vertex_const_iter Vertex_const_iter;
  typedef typename Mesh_properties::Face_list Face_list;
  typedef typename Mesh_properties::Face_iter Face_iter;
  typedef typename Mesh_properties::Face_const_iter Face_const_iter;
  // sample links
  typedef typename Mesh_properties::Link Link;
  typedef typename Mesh_properties::Link_list Link_list;
  typedef typename Mesh_properties::Link_list_iter Link_list_iter;
  typedef typename Mesh_properties::Link_list_const_iter Link_list_const_iter;
  typedef typename Mesh_properties::Link_iter_list Link_iter_list;
  typedef typename Mesh_properties::Link_iter_list_iter Link_iter_list_iter;
  typedef typename Mesh_properties::Link_iter_list_const_iter
                                    Link_iter_list_const_iter;
  typedef typename Mesh_properties::Link_pointer_list Link_pointer_list;
  typedef typename Mesh_properties::Link_pointer_iter Link_pointer_iter;
  typedef typename Mesh_properties::Link_pointer_const_iter
                                    Link_pointer_const_iter;
  // Point list
  typedef typename Mesh_properties::Point_list Point_list;
  typedef typename Mesh_properties::Point_iter Point_iter;
  typedef typename Mesh_properties::Point_const_iter Point_const_iter;
  // Color list
  typedef typename Mesh_properties::Color_list Color_list;
  typedef typename Mesh_properties::Color_iter Color_iter;
  typedef typename Mesh_properties::Color_const_iter Color_const_iter;
  // AABB tree
  typedef typename Mesh_properties::Face_tree Face_tree;
  // Dynamic priority queues
  typedef typename Mesh_properties::DPQueue_halfedge_long
                                    DPQueue_halfedge_long;
  typedef typename Mesh_properties::Halfedge_long Halfedge_long;
  typedef typename Mesh_properties::DPQueue_halfedge_short
                                    DPQueue_halfedge_short;
  typedef typename Mesh_properties::Halfedge_short Halfedge_short;
  typedef typename Mesh_properties::DPQueue_vertex_long DPQueue_vertex_long;
  typedef typename Mesh_properties::Vertex_long Vertex_long;
  typedef typename Mesh_properties::DPQueue_vertex_short DPQueue_vertex_short;
  typedef typename Mesh_properties::Vertex_short Vertex_short;
  typedef typename Mesh_properties::DPQueue_face_long DPQueue_face_long;
  typedef typename Mesh_properties::Face_long Face_long;
  typedef typename Mesh_properties::DPQueue_face_short DPQueue_face_short;
  typedef typename Mesh_properties::Face_short Face_short;
  // Visit list and iterator
  typedef typename std::list<std::pair<Point, FT>> Visit_list;
  typedef typename std::list<std::pair<Point, FT>>::iterator Visit_iter;

 public:
  // 1) life cycles
  Minangle_remesher() {
    // general paramters
    np_.max_error_threshold = 0.2;
    np_.min_angle_threshold = 30.0;
    np_.max_mesh_complexity = 100000000;
    np_.smooth_angle_delta = 0.1;
    np_.apply_edge_flip = true;
    np_.edge_flip_strategy = EdgeFlipStrategy::k_improve_angle;
    np_.flip_after_split_and_collapse = true;
    np_.relocate_after_local_operations = true;
    np_.relocate_strategy = RelocateStrategy::k_cvt_barycenter;
    np_.keep_vertex_in_one_ring = false;
    np_.use_local_aabb_tree = true;
    np_.collapsed_list_size = 10;
    np_.decrease_max_errors = true;
    np_.verbose_progress = true;
    np_.apply_initial_mesh_simplification = true;
    np_.apply_final_vertex_relocation = true;
    // sample parameters
    np_.samples_per_face_in = 10;
    np_.samples_per_face_out = 10;
    np_.max_samples_per_area = 10000;
    np_.min_samples_per_triangle = 1;
    np_.bvd_iteration_count = 1;
    np_.sample_number_strategy = SampleNumberStrategy::k_fixed;
    np_.sample_strategy = SampleStrategy::k_adaptive;
    np_.use_stratified_sampling = false;
    // feature parameters
    np_.sum_theta = 1.0;
    np_.sum_delta = 0.5;
    np_.dihedral_theta = 1.0;
    np_.dihedral_delta = 0.5;
    np_.feature_difference_delta = 0.15;
    np_.feature_control_delta = 0.5;
    np_.inherit_element_types = false;
    np_.use_feature_intensity_weights = false;
    // vertex optimization parameters
    np_.vertex_optimize_count = 2;
    np_.vertex_optimize_ratio = 0.9;
    np_.stencil_ring_size = 1;
    np_.optimize_strategy = OptimizeStrategy::k_approximation;
    np_.face_optimize_type = OptimizeType::k_both;
    np_.edge_optimize_type = OptimizeType::k_both;
    np_.vertex_optimize_type = OptimizeType::k_both;
    np_.optimize_after_local_operations = true;

    input_ = NULL;
    remesh_ = NULL;
    input_bbox = Bbox(DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX,
                      DOUBLE_MIN, DOUBLE_MIN, DOUBLE_MIN);
    links_initialized_ = false;
    input_aabb_tree_constructed_ = false;
  }

  explicit Minangle_remesher(const NamedParameters &np) {
    np_ = np;

    input_ = NULL;
    remesh_ = NULL;
    input_bbox = Bbox(DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX,
                      DOUBLE_MIN, DOUBLE_MIN, DOUBLE_MIN);
    links_initialized_ = false;
    input_aabb_tree_constructed_ = false;
  }

  virtual ~Minangle_remesher() {
    if (input_ != NULL) {
      delete input_;
      input_ = NULL;
    }
    if (remesh_ != NULL) {
      delete remesh_;
      remesh_ = NULL;
    }
  }

  // 2) parameter access functions
  // 2.1) general parameters
  FT get_max_error_threshold() const { return np_.max_error_threshold; }
  void set_max_error_threshold(FT value) { np_.max_error_threshold = value; }
  FT get_min_angle_threshold() const { return np_.min_angle_threshold; }
  void set_min_angle_threshold(FT value) { np_.min_angle_threshold = value; }
  int get_max_mesh_complexity() const { return np_.max_mesh_complexity; }
  void set_max_mesh_complexity(int value) { np_.max_mesh_complexity = value; }
  FT get_smooth_angle_delta() const { return np_.smooth_angle_delta; }
  void set_smooth_angle_delta(FT value) { np_.smooth_angle_delta = value; }
  bool get_apply_edge_flip() const { return np_.apply_edge_flip; }
  void set_apply_edge_flip(bool value) { np_.apply_edge_flip = value; }
  EdgeFlipStrategy get_edge_flip_strategy() const
      { return np_.edge_flip_strategy; }
  void set_edge_flip_strategy(EdgeFlipStrategy value)
      { np_.edge_flip_strategy = value; }
  bool get_flip_after_split_and_collapse() const
      { return np_.flip_after_split_and_collapse; }
  void set_flip_after_split_and_collapse(bool value)
      { np_.flip_after_split_and_collapse = value; }
  bool get_relocate_after_local_operations() const
      { return np_.relocate_after_local_operations; }
  void set_relocate_after_local_operations(bool value)
      { np_.relocate_after_local_operations = value; }
  RelocateStrategy get_relocate_strategy() const
      { return np_.relocate_strategy; }
  void set_relocate_strategy(RelocateStrategy value)
      { np_.relocate_strategy = value; }
  bool get_keep_vertex_in_one_ring() const
      { return np_.keep_vertex_in_one_ring; }
  void set_keep_vertex_in_one_ring(bool value)
      { np_.keep_vertex_in_one_ring = value; }
  bool get_use_local_aabb_tree() const { return np_.use_local_aabb_tree; }
  void set_use_local_aabb_tree(bool value) { np_.use_local_aabb_tree = value; }
  int get_collapsed_list_size() const { return np_.collapsed_list_size; }
  void set_collapsed_list_size(int value) { np_.collapsed_list_size = value; }
  bool get_decrease_max_errors() const { return np_.decrease_max_errors; }
  void set_decrease_max_errors(bool value) { np_.decrease_max_errors = value; }
  bool get_verbose_progress() const { return np_.verbose_progress; }
  void set_verbose_progress(bool value) { np_.verbose_progress = value; }
  bool get_apply_initial_mesh_simplification() const
      { return np_.apply_initial_mesh_simplification; }
  void set_apply_initial_mesh_simplification(bool value)
      { np_.apply_initial_mesh_simplification = value; }
  bool get_apply_final_vertex_relocation() const
      { return np_.apply_final_vertex_relocation; }
  void set_apply_final_vertex_relocation(bool value)
      { np_.apply_final_vertex_relocation = value; }
  // 2.2) sample parameters
  int get_samples_per_face_in() const { return np_.samples_per_face_in; }
  void set_samples_per_face_in(int value) { np_.samples_per_face_in = value; }
  int get_samples_per_face_out() const { return np_.samples_per_face_out; }
  void set_samples_per_face_out(int value)
      { np_.samples_per_face_out = value; }
  int get_max_samples_per_area() const { return np_.max_samples_per_area; }
  void set_max_samples_per_area(int value)
      { np_.max_samples_per_area = value; }
  int get_min_samples_per_triangle() const
      { return np_.min_samples_per_triangle; }
  void set_min_samples_per_triangle(int value)
      { np_.min_samples_per_triangle = value; }
  int get_bvd_iteration_count() const { return np_.bvd_iteration_count; }
  void set_bvd_iteration_count(int value) { np_.bvd_iteration_count = value; }
  SampleNumberStrategy get_sample_number_strategy() const
      { return np_.sample_number_strategy; }
  void set_sample_number_strategy(SampleNumberStrategy value)
      { np_.sample_number_strategy = value; }
  SampleStrategy get_sample_strategy() const { return np_.sample_strategy; }
  void set_sample_strategy(SampleStrategy value)
      { np_.sample_strategy = value; }
  bool get_use_stratified_sampling() const
      { return np_.use_stratified_sampling; }
  void set_use_stratified_sampling(bool value)
      { np_.use_stratified_sampling = value; }
  // 2.3) feature intensity parameters
  FT get_sum_theta() const { return np_.sum_theta; }
  void set_sum_theta(FT value) { np_.sum_theta = value; }
  FT get_sum_delta() const { return np_.sum_delta; }
  void set_sum_delta(FT value) { np_.sum_delta = value; }
  FT get_dihedral_theta() const { return np_.dihedral_theta; }
  void set_dihedral_theta(FT value) { np_.dihedral_theta = value; }
  FT get_dihedral_delta() const { return np_.dihedral_delta; }
  void set_dihedral_delta(FT value) { np_.dihedral_delta = value; }
  FT get_feature_difference_delta() const
      { return np_.feature_difference_delta; }
  void set_feature_difference_delta(FT value)
      { np_.feature_difference_delta = value; }
  FT get_feature_control_delta() const { return np_.feature_control_delta; }
  void set_feature_control_delta(FT value)
      { np_.feature_control_delta = value; }
  bool get_inherit_element_types() const { return np_.inherit_element_types; }
  void set_inherit_element_types(bool value)
      { np_.inherit_element_types = value; }
  bool get_use_feature_intensity_weights() const
      { return np_.use_feature_intensity_weights; }
  void set_use_feature_intensity_weights(bool value)
      { np_.use_feature_intensity_weights = value; }
  // 2.4) vertex optimization parameters
  int get_vertex_optimize_count() const { return np_.vertex_optimize_count; }
  void set_vertex_optimize_count(int value)
      { np_.vertex_optimize_count = value; }
  FT get_vertex_optimize_ratio() const { return np_.vertex_optimize_ratio; }
  void set_vertex_optimize_ratio(FT value)
      { np_.vertex_optimize_ratio = value; }
  int get_stencil_ring_size() const { return np_.stencil_ring_size; }
  void set_stencil_ring_size(int value) { np_.stencil_ring_size = value; }
  OptimizeStrategy get_optimize_strategy() const
      { return np_.optimize_strategy; }
  void set_optimize_strategy(OptimizeStrategy value)
      { np_.optimize_strategy = value; }
  OptimizeType get_face_optimize_type() const
      { return np_.face_optimize_type; }
  void set_face_optimize_type(OptimizeType value)
      { np_.face_optimize_type = value; }
  OptimizeType get_edge_optimize_type() const { return np_.edge_optimize_type; }
  void set_edge_optimize_type(OptimizeType value)
      { np_.edge_optimize_type = value; }
  OptimizeType get_vertex_optimize_type() const
      { return np_.vertex_optimize_type; }
  void set_vertex_optimize_type(OptimizeType value)
      { np_.vertex_optimize_type = value; }
  bool get_optimize_after_local_operations() const
      { return np_.optimize_after_local_operations; }
  void set_optimize_after_local_operations(bool value)
      { np_.optimize_after_local_operations = value; }

  // 3) member data access
  Bbox get_input_bbox() const { return input_bbox; }
  const NamedParameters &get_named_parameters() const { return np_; }
  bool get_links_initialized() const { return links_initialized_; }
  void set_input(Mesh *input, bool verbose_progress) {
    // step 1: set the input
    delete_input();
    input_ = new Mesh_properties(input);
    input_bbox = input_->calculate_bounding_box();
    // step 2: calculate the properties
    calculate_normals(true, verbose_progress);
    // step 3: update feature intensities and clear links
    calculate_feature_intensities(true, false, verbose_progress);
    // step 4: update status
    input_aabb_tree_constructed_ = false;
  }
  Mesh_properties* get_input() { return input_; }
  const Mesh_properties* get_input() const { return input_; }
  void set_remesh(Mesh *remesh, bool verbose_progress) {
    // step 1: set the remesh
    delete_remesh();
    remesh_ = new Mesh_properties(remesh);
    // step 2: calculate the properties
    calculate_normals(false, verbose_progress);
    // step 3: update feature intensities and clear links
    calculate_feature_intensities(false, true, verbose_progress);
  }
  Mesh_properties* get_remesh() { return remesh_; }
  const Mesh_properties* get_remesh() const { return remesh_; }
  void delete_input() {
    if (input_ != NULL) {
      delete input_;
    }
    input_ = NULL;
  }
  void delete_remesh() {
    if (remesh_ != NULL) {
      delete remesh_;
    }
    remesh_ = NULL;
  }
  void save_remesh_as(const std::string &file_name) const {
    size_t pos = file_name.find_last_of('.');
    if (pos == std::string::npos) {
      std::cout << "Invalid file name" << std::endl;
      return;
    }
    std::string extension = file_name.substr(pos);
    std::transform(extension.begin(), extension.end(),
      extension.begin(), std::tolower);
    if (extension == ".off") {
      save_as_off(file_name);
    } else {
      std::cout << "Invalid file type" << std::endl;
    }
  }

  // 4) mesh properties
  void input_properties() const {
    std::cout << std::endl;
    std::cout << yellow << "INPUT PROPERTIES" << white << std::endl;
    input_->trace_properties();
  }
  bool remesh_properties() {
    bool status_before = links_initialized_, status_after = links_initialized_;
    if (remesh_ != NULL) {
      if (!links_initialized_) {
        generate_samples_and_links();
        status_after = true;
      }
      std::cout << std::endl;
      std::cout << yellow << "REMESH PROPERTIES" << white << std::endl;
      remesh_->trace_properties();
      FT diagonal_length = input_->calculate_diagonal_length();
      remesh_->trace_additional_properties(diagonal_length);
    }
    return status_before != status_after;
  }

  // 5) max error
  FT get_max_error_threshold_value() const {
    FT diagonal = std::sqrt(std::pow(
      input_bbox.xmax() - input_bbox.xmin(), 2) +
      std::pow(input_bbox.ymax() - input_bbox.ymin(), 2) +
      std::pow(input_bbox.zmax() - input_bbox.zmin(), 2));
    return diagonal * np_.max_error_threshold * 0.01;
  }

  // 6) feature intentisities
  void calculate_feature_intensities(bool update_input, bool update_remesh,
    bool verbose_progress) {
    if (update_input && input_ != NULL) {
      if (verbose_progress) {
        std::cout << "Computing input feature intensities...";
      }
      input_->calculate_feature_intensities(np_);
      if (verbose_progress) {
        std::cout << "Done" << std::endl;
      }
    }
    if (update_remesh && remesh_ != NULL) {
      if (verbose_progress) {
        std::cout << "Computing remesh feature intensities...";
      }
      remesh_->calculate_feature_intensities(np_);
      if (verbose_progress) {
        std::cout << "Done" << std::endl;
      }
    }
    links_initialized_ = false;
  }

  // 7) sample and links
  void generate_samples_and_links() {
    if (input_ == NULL || remesh_ == NULL || links_initialized_) {
      return;
    }
    // step 1: generate AABB trees if necessary
    std::cout << std::endl;
    if (!input_aabb_tree_constructed_) {
      build_face_tree(true, &input_face_tree_);
      input_aabb_tree_constructed_ = true;
    }
    build_face_tree(false, &remesh_face_tree_);
    // step 2: clear all the links in input_ and remesh_
    clear_links();
    // step 3: generate the out links
    int samples_per_face = np_.samples_per_face_out;
    if (np_.sample_number_strategy == SampleNumberStrategy::k_variable) {
      FT value = static_cast<FT>(np_.samples_per_face_out);
      value *= input_->get_mesh().number_of_faces();
      value /= remesh_->get_mesh().number_of_faces();
      samples_per_face = static_cast<int>(value);
    }
    remesh_->generate_out_links(input_face_tree_, samples_per_face,
      INITIAL_BVD_COUNT, NULL, np_);
    // step 4: generate the in links
    input_->generate_out_links(remesh_face_tree_, np_.samples_per_face_in,
      INITIAL_BVD_COUNT, remesh_, np_);
    // step 5: compute the max_squared_errors
    remesh_->calculate_max_squared_errors();
    links_initialized_ = true;
  }

  // 8) polyhedron manipulations
  int eliminate_input_degenerated_faces() const {
    return input_->eliminate_degenerated_faces();
  }

  int split_input_long_edges() const {
    return input_->split_long_edges();
  }

  void minangle_remeshing() {
    if (!links_initialized_) {
      generate_samples_and_links();
    }
    CGAL::Timer timer;
    timer.start();
    std::cout << std::endl << "Min angle remeshing..." << std::endl;
    if (np_.apply_initial_mesh_simplification) {
      std::cout << std::endl;
      initial_mesh_simplification();
    }
    std::cout << std::endl;
    maximize_minimal_angle();
    if (np_.apply_final_vertex_relocation) {
      std::cout << std::endl;
      final_vertex_relocation();
    }
    std::cout << std::endl;
    std::cout << "Done, (total time is " << timer.time() << " s)" << std::endl;
  }

  void initial_mesh_simplification() {
    /* for isotropic purpose, fill priority queue with collapsible edges.
    The priority has two choices:
    1) the error before collapse: E_{before}
    2) edge_length * opposite_angle (current implementation) */
    if (!links_initialized_) {
      generate_samples_and_links();
    }
    FT max_error_threshold_value = get_max_error_threshold_value();
    CGAL::Timer timer;
    timer.start();
    std::cout << std::endl << "Initial mesh simplification..." << std::endl;
    std::cout << "(max error threshold value = "
      << max_error_threshold_value << ")" << std::endl;
    DPQueue_halfedge_long large_error_queue;
    DPQueue_halfedge_short collapse_candidate_queue;
    remesh_->fill_collapse_candidate_edges(max_error_threshold_value,
      &large_error_queue, &collapse_candidate_queue, np_);
    unsigned int index = 0, nb_operations = 0;
    halfedge_descriptor max_error_halfedge;
    FT max_error = 0.0;
    while (remesh_->size_of_vertices() < np_.max_mesh_complexity &&
      (!large_error_queue.empty() || !collapse_candidate_queue.empty())) {
      while (remesh_->size_of_vertices() < np_.max_mesh_complexity &&
        !large_error_queue.empty()) {
        DPQueue_halfedge_long::right_map::iterator eit =
          large_error_queue.right.begin();
        max_error = CGAL::sqrt(eit->first);
        max_error_halfedge = eit->second;
        large_error_queue.right.erase(eit);
        if (np_.verbose_progress) {
          std::cout << ++index << ": error queue size = "
            << large_error_queue.size() << " ";
        }
        greedy_reduce_error(max_error_threshold_value, max_error,
          np_.verbose_progress, true, &large_error_queue,
          &collapse_candidate_queue, max_error_halfedge);
      }
      if (!collapse_candidate_queue.empty()) {
        if (np_.verbose_progress) {
          std::cout << ++index << ": collapse queue size = "
            << collapse_candidate_queue.size() << " ";
        }
        // step 1: get the top halfedge that might be collapsed
        DPQueue_halfedge_short::right_map::iterator eit =
          collapse_candidate_queue.right.begin();
        halfedge_descriptor hd = eit->second;
        collapse_candidate_queue.right.erase(eit);
        collapse_candidate_queue.left.erase(remesh_->get_opposite(hd));
        // step 2: try to collapse with the constraints of max_error
        vertex_descriptor vd = collapse_applied(max_error_threshold_value,
          -1.0, true, NULL, &large_error_queue,
          &collapse_candidate_queue, hd);
        // step 3: update the adjacent halfedges
        if (vd != remesh_->get_null_vertex()) {
          ++nb_operations;
          if (np_.verbose_progress) {
            std::cout << "1 edge collapsed";
          }
        }
        if (np_.verbose_progress) {
          std::cout << std::endl;
        }
      }
    }
    std::cout << "Done (" << nb_operations << " local operations applied, "
      << timer.time() << " s)" << std::endl;
  }

  void split_local_longest_edge() {
    if (!links_initialized_) {
      generate_samples_and_links();
    }
    FT max_error = 0.0, min_radian = CGAL_PI;
    halfedge_descriptor max_error_halfedge, min_radian_halfedge;
    FT max_error_threshold_value = get_max_error_threshold_value();
    max_error_halfedge = remesh_->calculate_maximal_error(&max_error);
    min_radian_halfedge = remesh_->calculate_minimal_radian(&min_radian);
    std::cout << std::endl << "Split local longest edge..." << std::endl;
    std::cout << "(max error threshold value = " << max_error_threshold_value
      << ", min angle threshold = " << np_.min_angle_threshold
      << " degree)" << std::endl;
    face_descriptor fd = remesh_->get_face(min_radian_halfedge);
    halfedge_descriptor longest_hd = remesh_->get_longest_halfedge(fd);
    longest_hd = remesh_->longest_side_propagation(longest_hd);
    vertex_descriptor vd = remesh_->split_edge(input_face_tree_,
      max_error_threshold_value, max_error, min_radian,
      false, NULL, NULL, longest_hd, np_);
    if (vd != remesh_->get_null_vertex()) {
      std::cout << "1 local longest edge splitted" << std::endl;
    } else {
      std::cout << "Error: no edge splitted" << std::endl;
    }
  }

  void increase_minimal_angle() {
    if (!links_initialized_) {
      generate_samples_and_links();
    }
    // step 1: try to decrease the max error if necessary
    FT max_error_threshold_value = get_max_error_threshold_value();
    std::cout << std::endl << "Increase minimal angle..." << std::endl;
    if (np_.decrease_max_errors) {
      FT max_error = 0.0;
      halfedge_descriptor max_error_halfedge;
      max_error_halfedge = remesh_->calculate_maximal_error(&max_error);
      if (max_error >= max_error_threshold_value) {
        std::cout << "(max error threshold = " << max_error_threshold_value
          << ", max error = " << max_error << ")" << std::endl;
        greedy_reduce_error(max_error_threshold_value, max_error, true, false,
          NULL, NULL, max_error_halfedge);
        return;
      }
    }
    // step 2: try to increase the min radian
    FT min_radian = CGAL_PI;
    halfedge_descriptor min_radian_halfedge;
    min_radian_halfedge = remesh_->calculate_minimal_radian(&min_radian);
    std::cout << "(min angle threshold = " << np_.min_angle_threshold
      << " degree, min angle = " << remesh_->to_angle(min_radian)
      << " degree)" << std::endl;
    greedy_improve_angle(max_error_threshold_value, min_radian, true, NULL,
      NULL, min_radian_halfedge);
  }

  void maximize_minimal_angle() {
    if (!links_initialized_) {
      generate_samples_and_links();
    }
    FT max_error_threshold_value = get_max_error_threshold_value();
    FT min_radian_threshold = remesh_->to_radian(np_.min_angle_threshold);
    CGAL::Timer timer;
    timer.start();
    std::cout << std::endl << "Greedy angles improvement..." << std::endl;
    std::cout << "(max error threshold value = " << max_error_threshold_value
      << ", min angle threshold = " << np_.min_angle_threshold
      << " degree, max mesh complexity = " << np_.max_mesh_complexity
      << ")" << std::endl;
    FT max_error = 0, min_radian = CGAL_PI;
    halfedge_descriptor max_error_halfedge, min_radian_halfedge;
    unsigned int nb_operations = 0;
    // Version 1: use dynamic priority queue
    DPQueue_halfedge_long large_error_queue;
    DPQueue_halfedge_short small_radian_queue;
    remesh_->fill_small_radian_edges(max_error_threshold_value,
      &large_error_queue, &small_radian_queue, np_);
    while (remesh_->size_of_vertices() < np_.max_mesh_complexity &&
      (!large_error_queue.empty() || !small_radian_queue.empty())) {
      while (remesh_->size_of_vertices() < np_.max_mesh_complexity &&
        !large_error_queue.empty()) {
        DPQueue_halfedge_long::right_map::iterator eit =
          large_error_queue.right.begin();
        max_error = CGAL::sqrt(eit->first);
        max_error_halfedge = eit->second;
        large_error_queue.right.erase(eit);
        if (np_.verbose_progress) {
          std::cout << ++nb_operations << ": max error = "
            << max_error << " ";
        }
        greedy_reduce_error(max_error_threshold_value, max_error,
          np_.verbose_progress, false, &large_error_queue,
          &small_radian_queue, max_error_halfedge);
      }
      if (remesh_->size_of_vertices() < np_.max_mesh_complexity &&
        !small_radian_queue.empty()) {
        DPQueue_halfedge_short::right_map::iterator eit =
          small_radian_queue.right.begin();
        min_radian = eit->first;
        min_radian_halfedge = eit->second;
        small_radian_queue.right.erase(eit);
        if (np_.verbose_progress) {
          std::cout << ++nb_operations << ": min angle = "
            << remesh_->to_angle(min_radian) << " ";
        }
        greedy_improve_angle(max_error_threshold_value, min_radian,
          np_.verbose_progress, &large_error_queue, &small_radian_queue,
          min_radian_halfedge);
      }
    }
    // Version 2: do not use dynamic priority queue
    /*if (decrease_max_errors_) {
    max_error_halfedge = remesh_->get_maximal_error(&max_error);
    min_radian_halfedge = remesh_->get_minimal_radian(&min_radian);
    while (remesh_->size_of_vertices() < np_.max_mesh_complexity &&
    (max_error >= max_error_threshold ||
    min_radian < min_radian_threshold)) {
    while (remesh_->size_of_vertices() < np_.max_mesh_complexity &&
    max_error >= max_error_threshold) {
    if (np_.verbose_progress) {
    std::cout << ++nb_operations << ": max error = "
    << max_error << " ";
    }
    greedy_reduce_error(input_face_tree, max_error_threshold,
    max_error, np_.verbose_progress, false, NULL, NULL,
    max_error_halfedge, m_pRemesh);
    max_error_halfedge = remesh_->get_maximal_error(&max_error);
    min_radian_halfedge = remesh_->get_minimal_radian(&min_radian);
    }
    if (remesh_->size_of_vertices() < np_.max_mesh_complexity &&
    min_radian < min_radian_threshold) {
    if (np_.verbose_progress) {
    std::cout << ++nb_operations << ": min angle = "
    << remesh_->to_angle(min_radian) << " ";
    }
    greedy_improve_angle(input_face_tree, max_error_threshold,
    min_radian, np_.verbose_progress, NULL, NULL,
    min_radian_halfedge, m_pRemesh);
    max_error_halfedge = remesh_->get_maximal_error(&max_error);
    min_radian_halfedge = remesh_->get_minimal_radian(&min_radian);
    }
    }
    } else {
    min_radian_halfedge = remesh_->get_minimal_radian(&min_radian);
    while (remesh_->size_of_vertices() < np_.max_mesh_complexity &&
    min_radian < min_radian_threshold) {
    if (np_.verbose_progress) {
    std::cout << ++nb_operations << ": min angle = "
    << remesh_->to_angle(min_radian) << " ";
    }
    greedy_improve_angle(input_face_tree, max_error_threshold,
    min_radian, np_.verbose_progress, NULL, NULL,
    min_radian_halfedge, m_pRemesh);
    min_radian_halfedge = remesh_->get_minimal_radian(&min_radian);
    }
    }*/
    std::cout << "Done (" << nb_operations << " local operations applied, "
      << timer.time() << " s)" << std::endl;
  }

  void final_vertex_relocation() {
    if (!links_initialized_) {
      generate_samples_and_links();
    }
    FT max_error_threshold_value = get_max_error_threshold_value();
    CGAL::Timer timer;
    timer.start();
    std::cout << std::endl << "Final vertex relocation..." << std::endl;
    std::cout << "(max error threshold value = " << max_error_threshold_value
      << ", min angle threshold = " << np_.min_angle_threshold
      << " degree)..." << std::endl;
    DPQueue_vertex_short relocate_candidate_queue;
    remesh_->fill_relocate_candidate_vertices(&relocate_candidate_queue);
    unsigned int index = 0, nb_relocate = 0;
    while (!relocate_candidate_queue.empty()) {
      if (np_.verbose_progress) {
        std::cout << ++index << ": relocate queue size = "
          << relocate_candidate_queue.size() << " ";
      }
      // step 1: get the top vertex that might be relocated
      DPQueue_vertex_short::right_map::iterator eit =
        relocate_candidate_queue.right.begin();
      FT min_radian = eit->first;
      vertex_descriptor vd = eit->second;
      relocate_candidate_queue.right.erase(eit);
      // step 2: try to relocate with the constrait of max_error and min_radian
      Point initial_point = remesh_->calculate_initial_point_for_relocate(
        input_face_tree_, vd, np_);
      bool relocated = remesh_->relocate_vertex(input_face_tree_,
        max_error_threshold_value, -1.0, min_radian, false, NULL, NULL,
        initial_point, vd, np_);
      if (relocated) {
        remesh_->update_relocate_candidate_vertices(
          vd, &relocate_candidate_queue);
        ++nb_relocate;
        if (np_.verbose_progress) {
          std::cout << "1 vertices relocated";
        }
      }
      if (np_.verbose_progress) {
        std::cout << std::endl;
      }
    }
    std::cout << "Done (" << nb_relocate << " vertices relocated, "
      << timer.time() << " s)" << std::endl;
  }

 private:
  // 1) normals
  void calculate_normals(bool is_input, bool verbose_progress) const {
    const std::string name = is_input ? "input" : "remesh";
    if (verbose_progress) {
      std::cout << "Computing " + name + " normals...";
    }
    if (is_input) {
      input_->calculate_normals();
    } else {
      remesh_->calculate_normals();
    }
    if (verbose_progress) {
      std::cout << "Done" << std::endl;
    }
  }

  // 2) trees
  void build_face_tree(bool is_input, Face_tree *face_tree) const {
    CGAL::Timer timer;
    timer.start();
    if (is_input) {
      std::cout << "Building input face AABB tree...";
      input_->build_face_tree(face_tree);
    } else {
      std::cout << "Building remesh face AABB tree...";
      remesh_->build_face_tree(face_tree);
    }
    std::cout << "done (" << timer.time() << " s)" << std::endl;
  }

  // 3) IO
  void save_as_off(const std::string &file_name) const {
    if (remesh_ == NULL) {
      std::cout << "Please set the remesh first" << std::endl;
      return;
    }
    std::ofstream ofs(file_name);
    if (!ofs) {
      std::cout << "Failed to create the file" << std::endl;
      return;
    }
    const Mesh &mesh = remesh_->get_mesh();
    bool ret = write_off(ofs, mesh);
    if (!ret) {
      std::cout << "Failed to save the remesh" << std::endl;
    }
    ofs.close();
  }

  // 4) sample and links
  void clear_links() {
    // step 1: clear the out links
    remesh_->clear_out_links();
    // step 2: clear the in links (out links from the perspective of m_pInput)
    remesh_->clear_in_link_iterators();
    input_->clear_out_links();
    // step 3: clear private data
    collapsed_list_.clear();
    collapsed_map_.clear();
  }

  // 5) manipulations
  void greedy_improve_angle(FT max_error_threshold_value, FT min_radian,
      bool verbose_progress, DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_radian_queue,
      halfedge_descriptor min_radian_halfedge) {
    // improve min_radian by the following:
    // 1) If collapse applies, collapse and return;
    // 2) If flip applies, flip and return;
    // 3) If relocate applies, relocate and return;
    // 4) Find the local_longest_hh. If flip applies, flip local_longest_hh;
    //                               Otherwwiese, split local_longest_hh.
    // step 1: try to collapse
    bool infinite_loop = false;
    vertex_descriptor vd = collapse_applied(max_error_threshold_value,
        min_radian, false, &infinite_loop, large_error_queue,
        small_radian_queue, min_radian_halfedge);
    if (vd != remesh_->get_null_vertex()) {
      if (verbose_progress) {
        std::cout << "1 edge collapsed" << std::endl;
      }
      return;
    }
    // if no infinite loop encountered, try flip and relocate
    if (!infinite_loop) {
      // step 2: try to flip
      if (np_.apply_edge_flip) {
        int nb_flip = remesh_->flip_applied(input_face_tree_,
            max_error_threshold_value, -1.0, min_radian, false,
            large_error_queue, small_radian_queue, min_radian_halfedge, np_);
        if (nb_flip > 0) {
          if (verbose_progress) {
            std::cout << nb_flip << " edges flipped" << std::endl;
          }
          return;
        }
      }
      // step 3: try to relocate
      int nb_relocate = remesh_->relocate_applied(input_face_tree_,
          max_error_threshold_value, -1.0, min_radian, false,
          large_error_queue, small_radian_queue, min_radian_halfedge, np_);
      if (nb_relocate > 0) {
        if (verbose_progress) {
          std::cout << nb_relocate << " vertices relocated" << std::endl;
        }
        return;
      }
    }
    // step 3: split for later improvement
    face_descriptor fh = remesh_->get_face(min_radian_halfedge);
    halfedge_descriptor longest_hh = remesh_->get_longest_halfedge(fh);
    longest_hh = remesh_->longest_side_propagation(longest_hh);
    if (np_.apply_edge_flip) {
      halfedge_descriptor hnew = remesh_->flip_edge(input_face_tree_,
          max_error_threshold_value, -1.0, min_radian, false,
          large_error_queue, small_radian_queue, longest_hh, np_);
      if (hnew != remesh_->get_null_halfedge()) {
        if (verbose_progress) {
          std::cout << "1 longest edge flipped" << std::endl;
        }
        return;
      }
    }
    remesh_->split_edge(input_face_tree_, max_error_threshold_value, -1.0,
        min_radian, false, large_error_queue, small_radian_queue,
        longest_hh, np_);
    if (verbose_progress) {
      std::cout << "1 edge splitted" << std::endl;
    }
  }

  void greedy_reduce_error(FT max_error_threshold, FT max_error,
      bool verbose_progress, bool reduce_complexity,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue,
      halfedge_descriptor max_error_halfedge) const {
    // 1. reduce max_error in the following order:
    // 1) If flip reduces the max_error, flip and return;
    // 2) If relocate reduces the max_error, relocate and return;
    // 3) Otherwise, split the max_error_halfedge.
    // 2. if improve_min_radian is true, small_value_queue is the
    //    small_radian_queue; otherwise, it is the collapse_candidate_queue.
    // step 1: try to flip
    if (np_.apply_edge_flip && remesh_->flip_edge(input_face_tree_, -1.0,
        max_error, -1.0, reduce_complexity, large_error_queue,
        small_value_queue, max_error_halfedge, np_) !=
        remesh_->get_null_halfedge()) {
      if (verbose_progress) {
        std::cout << "1 edge flipped" << std::endl;
      }
      return;
    }
    // step 2: try to relocate
    int nb_relocate = remesh_->relocate_applied(input_face_tree_, -1.0,
        max_error, -1.0, reduce_complexity, large_error_queue,
        small_value_queue, max_error_halfedge, np_);
    if (nb_relocate > 0) {
      if (verbose_progress) {
        std::cout << nb_relocate << " edges relocated" << std::endl;
      }
      return;
    }
    // step 3: try to split
    remesh_->split_edge(input_face_tree_, -1.0, max_error, -1.0,
        reduce_complexity, large_error_queue, small_value_queue,
        max_error_halfedge, np_);
    if (verbose_progress) {
      std::cout << "1 edge splitted" << std::endl;
    }
  }

  bool caused_infinite_loop(halfedge_descriptor hd) {
    std::map<Point, std::map<FT, Visit_iter>, Point_Comp>::iterator it1;
    Visit_iter it2;
    Point point = remesh_->get_point(remesh_->get_opposite_vertex(hd));
    FT sl = remesh_->squared_length(hd);
    sl = to_approximation(sl);
    bool found = false;
    it1 = collapsed_map_.find(point);
    if (it1 != collapsed_map_.end()) {
      std::map<FT, Visit_iter>::iterator it = it1->second.find(sl);
      if (it != it1->second.end()) {
        found = true;
        it2 = it->second;
      }
    }
    if (found) {
      if (np_.verbose_progress) {
        std::cout << "Point(" << point << ") with length "
          << sl << ": collapse denied.";
      }
      // option: we may delete it instead of splicing front
      collapsed_list_.splice(collapsed_list_.begin(),
        collapsed_list_, it2);
      return true;
    } else {
      collapsed_list_.push_front(std::pair<Point, FT>(point, sl));
      collapsed_map_[point][sl] = collapsed_list_.begin();
      while (collapsed_list_.size() > np_.collapsed_list_size) {
        const Point &p = collapsed_list_.back().first;
        FT sl = collapsed_list_.back().second;
        it1 = collapsed_map_.find(p);
        if (it1 != collapsed_map_.end()) {
          it1->second.erase(sl);
          if (it1->second.size() == 0) {
            collapsed_map_.erase(it1);
          }
        }
        collapsed_list_.pop_back();
      }
      return false;
    }
  }

  // 6) collapse
  vertex_descriptor collapse_applied(FT max_error_threshold_value,
      FT min_radian, bool reduce_complexity, bool *infinite_loop,
      DPQueue_halfedge_long *large_error_queue,
      DPQueue_halfedge_short *small_value_queue, halfedge_descriptor hd) {
    /* if min_radian > 0, we improve min_radian;
    if infinite_loop is not NULL, we check the infinite loop case;
    if improve_min_radian, we improve the min radian; otherwise,
    we collapse to reduce the mesh complexity */
    // step 1: topology constraints check
    if (!remesh_->is_collapsable(hd)) {
      return remesh_->get_null_vertex();
    }
    // step 2: geometry constraints check
    Halfedge_list halfedges;  // use the halfedges to represent faces
    bool is_ring = remesh_->predict_faces_after_collapse(hd, &halfedges);
    Point new_point = remesh_->calculate_initial_point_for_collapse(hd, np_);
    if (np_.keep_vertex_in_one_ring &&
        remesh_->collapse_would_cause_wrinkle(halfedges, new_point, hd)) {
      return remesh_->get_null_vertex();
    }
    // step 3: backup the original local links
    std::set<face_descriptor> one_ring_faces, extended_faces;
    remesh_->collect_one_ring_faces_incident_to_edge(hd, &one_ring_faces);
    remesh_->extend_faces(one_ring_faces, np_.stencil_ring_size,
                          &extended_faces);
    Link_iter_list face_in_links, edge_in_links;
    Link_pointer_list vertex_in_links;
    Point_list face_in_end_points, edge_in_end_points, vertex_in_end_points;
    remesh_->backup_local_in_links(extended_faces, &face_in_links,
        &face_in_end_points, &edge_in_links, &edge_in_end_points,
        &vertex_in_links, &vertex_in_end_points);
    // step 4: simulate the edge collapse
    FT error = DOUBLE_MAX, radian = 0.0;
    simulate_edge_collapse(one_ring_faces, extended_faces, halfedges, hd,
        is_ring, face_in_links, edge_in_links, vertex_in_links, &error,
        &radian, &new_point);
    remesh_->restore_local_in_links(face_in_end_points, face_in_links,
        edge_in_end_points, edge_in_links, vertex_in_end_points,
        vertex_in_links);
    // step 5: fidelity constraints check (max_error)
    if (error >= max_error_threshold_value) {
      return remesh_->get_null_vertex();
    }
    // step 6: quality constraints check (min_radian) if necessary
    if (min_radian > 0 && radian < min_radian) {
      return remesh_->get_null_vertex();
    }
    // step 7: infinite loops case check if necessary
    if (infinite_loop != NULL) {
      *infinite_loop = caused_infinite_loop(hd);
      if (*infinite_loop) {
        return remesh_->get_null_vertex();
      }
    }
    // step 8: collapse the edge authentically
    vertex_descriptor vh = remesh_->collapse_edge(input_face_tree_,
      max_error_threshold_value, min_radian, reduce_complexity,
      large_error_queue, small_value_queue, face_in_links, edge_in_links,
      vertex_in_links, hd, new_point, np_);
    return vh;
  }

  void simulate_edge_collapse(const std::set<face_descriptor> &one_ring_faces,
      const std::set<face_descriptor> &extended_faces,
      const Halfedge_list &halfedges, halfedge_descriptor hh, bool is_ring,
      const Link_iter_list &face_in_links, const Link_iter_list &edge_in_links,
      const Link_pointer_list &vertex_in_links, FT *error, FT *radian,
      Point *new_point) const {
    // step 1: construct the local_mesh
    Mesh local_mesh;
    vertex_descriptor local_vd = remesh_->construct_local_mesh(one_ring_faces,
        extended_faces, halfedges, *new_point, is_ring, &local_mesh);
    Mesh_properties local_mp(&local_mesh);
    local_mp.calculate_feature_intensities(np_);
    // step 2: get the in_link_faces (for function compatability)
    std::set<face_descriptor> in_link_faces;
    local_mp.collect_all_faces(&in_link_faces);
    local_mp.generate_local_links(input_face_tree_, true, face_in_links,
        edge_in_links, vertex_in_links, local_vd, in_link_faces, np_);
    // step 3: optimize the vertex position if necessary
    if (np_.optimize_after_local_operations) {
      local_mp.optimize_vertex_position(input_face_tree_, face_in_links,
          edge_in_links, vertex_in_links, local_vd, in_link_faces, np_);
    }
    // step 4: update the max_errors for faces
    local_mp.calculate_max_squared_errors(&in_link_faces);
    // step 5: calculate the error, radian and new_point
    local_mp.calculate_local_maximal_error(in_link_faces, error);
    *radian = local_mp.calculate_minimal_radian_around_vertex(local_vd);
    *new_point = local_mesh.point(local_vd);
  }

  // 7) utilities
  inline FT to_approximation(FT value) const {
    FT precison = MAX_VALUE;
    int temp_value = value * precison;
    return temp_value / precison;
  }

 private:
  // 1) parameters
  NamedParameters np_;

  // 2) the collapse operator
  Visit_list collapsed_list_;
  std::map<Point, std::map<FT, Visit_iter>, Point_Comp> collapsed_map_;

  // 3) member data and properties
  Mesh_properties *input_, *remesh_;
  Face_tree input_face_tree_, remesh_face_tree_;
  Bbox input_bbox;

  // 4) status data
  bool links_initialized_;
  bool input_aabb_tree_constructed_;

  // 5) const data
  int const INITIAL_BVD_COUNT = 5;
};

}  // namespace internal
}  // namespace Polygon_mesh_processing
}  // namespace CGAL

#endif  // SRC_INTERNAL_MINANGLE_REMESHING_MINANGLE_REMESH_IMPL_H_
