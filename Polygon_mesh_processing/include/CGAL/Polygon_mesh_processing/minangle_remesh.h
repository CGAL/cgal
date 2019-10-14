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

#ifndef SRC_MINANGLE_REMESH_H_
#define SRC_MINANGLE_REMESH_H_

#include <string>
#include "internal\Minangle_remeshing\minangle_remesh_impl.h"

namespace CGAL {
namespace Polygon_mesh_processing {

template<class Kernel>
class Minangle_remesh {
 public:
  // typedef typenames
  typedef typename internal::Minangle_remesher<Kernel> Minangle_remesher;
  typedef typename Minangle_remesher::FT FT;
  typedef typename Minangle_remesher::Mesh Mesh;

  // life cycle
  Minangle_remesh() {
    remesher_ = new Minangle_remesher();
  }

  Minangle_remesh(FT max_error_threshold,     // general parameters
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
      bool verbose_progress,
      bool apply_initial_mesh_simplification,
      bool apply_final_vertex_relocation,
      int samples_per_face_in,                // sample parameters
      int samples_per_face_out,
      int max_samples_per_area,
      int min_samples_per_triangle,
      int bvd_iteration_count,
      SampleNumberStrategy sample_number_strategy,
      SampleStrategy sample_strategy,
      bool use_stratified_sampling,
      FT sum_theta,                           // feature intensity parameters
      FT sum_delta,
      FT dihedral_theta,
      FT dihedral_delta,
      FT feature_difference_delta,
      FT feature_control_delta,
      bool inherit_element_types,
      bool use_feature_intensity_weights,
      int vertex_optimize_count,              // vertex optimization parameters
      FT vertex_optimize_ratio,
      int stencil_ring_size,
      OptimizeStrategy optimize_strategy,
      OptimizeType face_optimize_type,
      OptimizeType edge_optimize_type,
      OptimizeType vertex_optimize_type,
      bool optimize_after_local_operations) {
    NamedParameters np;
    // general parameters
    np.max_error_threshold = max_error_threshold;
    np.min_angle_threshold = min_angle_threshold;
    np.max_mesh_complexity = max_mesh_complexity;
    np.smooth_angle_delta = smooth_angle_delta;
    np.apply_edge_flip = apply_edge_flip;
    np.edge_flip_strategy = edge_flip_strategy;
    np.flip_after_split_and_collapse = flip_after_split_and_collapse;
    np.relocate_after_local_operations = relocate_after_local_operations;
    np.relocate_strategy = relocate_strategy;
    np.keep_vertex_in_one_ring = keep_vertex_in_one_ring;
    np.use_local_aabb_tree = use_local_aabb_tree;
    np.collapsed_list_size = collapse_list_size;
    np.decrease_max_errors = decrease_max_errors;
    np.verbose_progress = verbose_progress;
    np.apply_initial_mesh_simplification = apply_initial_mesh_simplification;
    np.apply_final_vertex_relocation = apply_final_vertex_relocation;
    // sample parameters
    np.samples_per_face_in = samples_per_face_in;
    np.samples_per_face_out = samples_per_face_out;
    np.max_samples_per_area = max_samples_per_area;
    np.min_samples_per_triangle = min_samples_per_triangle;
    np.bvd_iteration_count = bvd_iteration_count;
    np.sample_number_strategy = sample_number_strategy;
    np.sample_strategy = sample_strategy;
    np.use_stratified_sampling = use_stratified_sampling;
    // feature intensity parameters
    np.sum_theta = sum_theta;
    np.sum_delta = sum_delta;
    np.dihedral_theta = dihedral_theta;
    np.dihedral_delta = dihedral_delta;
    np.feature_difference_delta = feature_difference_delta;
    np.feature_control_delta = feature_control_delta;
    np.inherit_element_types = inherit_element_types;
    np.use_feature_intensity_weights = use_feature_intensity_weights;
    // vertex optimization parameters
    np.vertex_optimize_count = vertex_optimize_count;
    np.vertex_optimize_ratio = vertex_optimize_ratio;
    np.stencil_ring_size = stencil_ring_size;
    np.optimize_strategy = optimize_strategy;
    np.face_optimize_type = face_optimize_type;
    np.edge_optimize_type = edge_optimize_type;
    np.vertex_optimize_type = vertex_optimize_type;
    np.optimize_after_local_operations = optimize_after_local_operations;

    remesher_ = new Minangle_remesher(np);
  }

  virtual ~Minangle_remesh() {
    free_minangle_remesher();
  }

  void free_minangle_remesher() {
    if (remesher_ != NULL) {
      delete remesher_;
      remesher_ = NULL;
    }
  }

  // public functions
  void set_input(Mesh *input, bool verbose_progress) const {
    remesher_->set_input(input, verbose_progress);
  }

  void set_remesh(Mesh *remesh, bool verbose_progress) const {
    remesher_->set_remesh(remesh, verbose_progress);
  }

  void save_remesh_as(const std::string &file_name) const {
    remesher_->save_remesh_as(file_name);
  }

  void delete_input() const { remesher_->delete_input(); }

  void delete_remesh() const { remesher_->delete_remesh(); }

  void generate_samples_and_links() const {
    remesher_->generate_samples_and_links();
  }

  void minangle_remeshing() const {
    remesher_->minangle_remeshing();
  }

  void initial_mesh_simplification() const {
    remesher_->initial_mesh_simplification();
  }

  void increase_minimal_angle() const {
    remesher_->increase_minimal_angle();
  }

  void maximize_minimal_angle() const {
    remesher_->maximize_minimal_angle();
  }

  void final_vertex_relocation() const {
    remesher_->final_vertex_relocation();
  }

  // access functions
  Minangle_remesher* get_remesher() { return remesher_; }

  const Minangle_remesher* get_remesher() const { return remesher_; }

 private:
  Minangle_remesher *remesher_;
};

}   // namespace Polygon_mesh_processing
}   // namespace CGAL

#endif  // SRC_MINANGLE_REMESH_H_
