// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Soichiro Yamazaki <soichiro19998@gmail.com>, Théo Bénard <benard320@gmail.com>
//
#ifndef HEXMESHING_FOR_LINEAR_CELL_COMPLEX_SEQUENTIAL_H
#define HEXMESHING_FOR_LINEAR_CELL_COMPLEX_SEQUENTIAL_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_prod_cons.h>
#include <CGAL/hexmeshing/Hexmeshing_grid.h>
#include <CGAL/hexmeshing/Hexmeshing_function_alias.h>
#include <CGAL/hexmeshing/Hexmeshing_function_generator.h>
#include <CGAL/hexmeshing/Hexmeshing_load_patterns.h>
#include <CGAL/hexmeshing/Hexmeshing_two_refinement_algorithm.h>
#include <CGAL/hexmeshing/Hexmeshing_post_processing.h>
#include <CGAL/query_replace/cmap_query_replace.h>
#include <CGAL/Hexmeshing_mesh_data_for_hexmeshing.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_traits.h>
#include <vector>



namespace CGAL {


/**
 * @brief Core data structure for hexahedral mesh generation and refinement
 *
 * This structure maintains the state of the hexahedral meshing process,
 * including the Linear Cell Complex (LCC), grid configuration, and various
 * markers used during mesh generation and refinement.
 */
class Hexmeshing_for_linear_cell_complex {
public:
  using Kernel = internal::Hexmeshing::Kernel;
  using LCC = internal::Hexmeshing::LCC;
  using Dart_handle = internal::Hexmeshing::Dart_handle;
  using DartInfo = internal::Hexmeshing::DartInfo;
  using size_type = internal::Hexmeshing::size_type;
  using PlaneCC = std::vector<Dart_handle>; // One dart per face connected components
  using PlaneSet = std::vector<PlaneCC>; // A set of planes
  /**
   * @brief Container for external pattern substitution resources used in hexahedral meshing
   *
   * This structure holds pattern substituters for both regular and partial templates.
   * These templates are used during the mesh refinement process to replace existing
   * mesh patterns with more refined ones.
   */
  struct ExternalRessources {
    Pattern_substituer<LCC> regular_templates;   ///< Pattern substituter for regular hexahedral templates
    Pattern_substituer<LCC> partial_templates;   ///< Pattern substituter for partial hexahedral templates
  };
  // Required initialization
  internal::Hexmeshing::Grid grid;                  ///< Grid configuration defining the mesh structure
  ExternalRessources* ext;    ///< Pointer to external resources for pattern substitution

  // Initialized by the algorithm
  LCC lcc;                    ///< Linear Cell Complex representing the mesh
  size_type identified_mark;  ///< Mark for identified cells
  size_type template_mark;    ///< Mark for template cells
  size_type propagation_face_mark;  ///< Mark for faces during propagation
  size_type debug;           ///< Debug marker
  size_type debug2;          ///< Secondary debug marker
  size_type three_template_node_mark;  ///< Mark only for ProcessData
  int level = 0;             ///< Current refinement level
  std::array<PlaneSet, 3> first_face_of_planes;  ///< First faces of each plane set (X, Y, Z)

  /// @brief Default constructor
  Hexmeshing_for_linear_cell_complex() {}

  /**
   * @brief Fixes invalid dart handles after refinement
   *
   * After refinement, some dart handles may become invalid. This function
   * iterates over all attributes to replace invalid handles with valid ones.
   * It ensures the consistency of the mesh structure after modifications.
   */
  void fix_dart_storage() {
    // Data that we'll use to know which plane/cc we are looking for, while iterating on attributes
    using DartRef = std::reference_wrapper<Dart_handle>;
    using CCIdToFace = std::unordered_map<size_t, DartRef>;
    using InvalidPlaneSet = std::unordered_map<size_t, CCIdToFace>;
    std::array<InvalidPlaneSet, 3> plane_non_valid_faces;

    const auto valid_face = [&](Dart_handle& face, char p, size_t cc_id){
      auto face_attr = lcc.attribute<2>(face);
      auto& face_info = face_attr->info();
      return lcc.is_dart_used(face) && face_attr != nullptr && face_attr->is_valid() && face_info.plane[p] == true && face_info.cc_id == cc_id;
    };

    bool should_fix = false;
    for (size_t p = 0; p < 3; p++){
      PlaneSet& plane_set = first_face_of_planes[p];
      for (size_t plane_id = 0; plane_id < plane_set.size(); plane_id++){
        PlaneCC& plane_cc = plane_set[plane_id];
        CCIdToFace non_valid_face;

        for (int cc_id = 0; cc_id < plane_cc.size(); cc_id++){
          Dart_handle& dart = plane_cc[cc_id];

          if (!valid_face(dart, p, cc_id)){
            non_valid_face.emplace(cc_id, std::ref(dart));
            should_fix = true;
          }
        }

        if (non_valid_face.size() > 0)
        plane_non_valid_faces[p][plane_id] = std::move(non_valid_face);
      }
    }

    if (!should_fix) return;

    // We don't care if a face has multiple planes, there is atleast one face for each cc that is not a merged face
    // Eitherway we cannot treat a face belonging to multiple planes, because the cc_id is undefined
    const auto get_plane_normal = [&](DartInfo::FaceAttrValue& face) -> std::optional<char>{
      std::bitset<3>& plane = face.plane;
      char sum = plane[0] + plane[1] + plane[2];
      if (sum > 1 or sum == 0) return {};

      for (char i = 0; i < 3; i++){
        if (plane[i]) return i;
      }

      CGAL_assertion_msg(false, "get_plane_normal: Unreachable code accessed");
      return {};
    };

    const auto fix_face_orientation = [&](Dart_handle face, char plane_normal) -> Dart_handle {
      auto n0 = lcc.attribute<0>(face)->point();
      auto n1 = lcc.attribute<0>(lcc.beta(face, 1))->point();
      auto n2 = lcc.attribute<0>(lcc.beta(face, 1, 1))->point();
      // auto n3 = lcc.attribute<0>(lcc.beta(face, 1, 1, 1));

      internal::Hexmeshing::Vector normal = cross_product(n1 - n0, n2 - n1);
      internal::Hexmeshing::Vector plane_up = [&]() -> internal::Hexmeshing::Vector {
        switch (plane_normal){
        default: return {1,0,0};
        case 1: return {0,1,0};
        case 2: return {0,0,1};
        }
      }();

      if (plane_up * normal < 0) {
        assert(!lcc.is_free<3>(face));
        return lcc.beta(face, 3);
      }

      return face;
    };

    // Invalidated handles are rare, we can afford scanning all 2 attributes
    auto& attributes = lcc.attributes<2>();
    for (auto it = attributes.begin(); it != attributes.end(); it++){
      auto& info = it->info();
      std::optional<char> normal = get_plane_normal(info);
      if (!normal) continue;

      auto& nvf_plane_set = plane_non_valid_faces[*normal];

      // Get plane
      auto plane_it = nvf_plane_set.find(info.plane_id);
      if (plane_it == nvf_plane_set.end()) continue;
      CCIdToFace& non_valid_faces = plane_it->second;

      // Get cc_id that was invalidated
      auto face_it = non_valid_faces.find(info.cc_id);
      if (face_it == non_valid_faces.end()) continue;

      // If both plane and cc_id match, reassign the handle
      DartRef handle_ptr = face_it->second;
      non_valid_faces.erase(face_it);

      handle_ptr.get() = fix_face_orientation(it->dart(), *normal);

      if (non_valid_faces.size() == 0){
        nvf_plane_set.erase(plane_it);
      }
    }

    auto all_valid = [&](){
      for (int p = 0; p < 3; p++){
        for (auto& non_valid : plane_non_valid_faces[p]){
          CGAL_assertion_msg(non_valid.second.size() == 0, "Hexmeshing_for_linear_cell_complex::fix_planes_sets, not all stored planes were fixed");
        }
      }

      for (int p = 0; p < 3; p++){
        auto& plane_set = first_face_of_planes[p];
        for (int plane_id = 0; plane_id < plane_set.size(); plane_id++){
          auto& plane_cc = plane_set[plane_id];
          for (int cc_id = 0; cc_id < plane_cc.size(); cc_id++){
            auto& face = plane_cc[cc_id];
            auto face_attr = lcc.attribute<2>(face);
            CGAL_assertion_msg(valid_face(face, p, cc_id), "Hexmeshing_for_linear_cell_complex::fix_planes_sets, face was not valid after fix");
          }
        }
      }
    };

    CGAL_assertion_code(all_valid());
  }

  virtual void reset_temp_vertex_ids() {}

  /**
   * @brief Initializes the mesh data with external resources and grid configuration
   * @param ext Pointer to external resources
   * @param grid Grid configuration for the mesh
   */
  void init(ExternalRessources* ext, internal::Hexmeshing::Grid grid) {
    this->ext = ext;
    this->grid = grid;
  }

  /**
   * @brief Main entry point for the two-refinement hexahedral mesh generation algorithm
   *
   * This function provides a high-level interface for generating hexahedral meshes
   * using the two-refinement algorithm. It orchestrates the complete process from
   * initialization to final mesh generation, including pattern loading, algorithm
   * execution, optional trimming, and post-processing with volume fraction analysis.
   *
   * The function performs the following operations:
   *
   * 1. **Resource Initialization**: Creates `Hexmeshing_for_linear_cell_complex` and `ExternalRessources`
   *    structures to hold the mesh data and pattern substitution resources
   *
   * 2. **Pattern Loading**: Calls `load_patterns` to load all necessary template
   *    patterns for both regular and partial template substitution
   *
   * 3. **Data Initialization**: Initializes the hexahedral meshing data with the
   *    provided grid configuration and external resources
   *
   * 4. **Algorithm Execution**: Calls `two_refinement_algorithm` to perform the
   *    complete refinement process with the specified number of levels
   *
   * 5. **Optional Trimming**: If trim is true, calls `trim_excedent_volumes` to remove
   *    excess volumes using the provided trimming function
   *
   * 6. **Post-Processing**: Applies volume fraction analysis and mesh smoothing:
   *    - Sets volume fractions based on the domain geometry
   *    - Performs volume fraction-based mesh optimization
   *    - Applies surface and volume smoothing for improved mesh quality
   *
   * The function returns a Linear Cell Complex containing the refined hexahedral mesh
   * with optimized geometry and volume fractions. The mesh quality is enhanced through
   * both topological refinement and geometric post-processing.
   *
   * @pre grid have cubic cells  (grid.size.x() = grid.size.y() = grid.size.z())
   * @param mesh Grid configuration defining the initial mesh structure and dimensions
   * @param nb_levels Number of refinement levels to perform (default: 1)
   * @param trim Whether to apply trimming to remove excess volumes after refinement
   *             (default: false). When true, trimmingFunction is used to determine
   *             which volumes to keep in the final mesh.
   */
  void two_refinement(
      Mesh_data_for_hexmeshing& mesh,
      int nb_levels = 1,
      bool trim = false)
  {
    using namespace internal::Hexmeshing;

    Tree* tree = mesh.get_tree_pointer();
    MarkingFunction cellIdentifier = is_volume_intersecting_poly(*tree);
    DecideInsideFunction decideFunc = is_inner_point(*tree);

    ExternalRessources res;

    load_patterns(res.regular_templates, res.partial_templates);
    init(&res, *mesh.get_grid_pointer());

    two_refinement_algorithm(*this, cellIdentifier, nb_levels);

    // assumes grid cells to be cubes
    post_processing(lcc, grid.size.x()/(1<<nb_levels), trim, cellIdentifier, decideFunc);
  }

  void two_refinement_without_post_processing(
      Mesh_data_for_hexmeshing& mesh,
      int nb_levels = 1)
  {
    using namespace internal::Hexmeshing;

    Tree* tree = mesh.get_tree_pointer();
    MarkingFunction cellIdentifier = is_volume_intersecting_poly(*tree);

    ExternalRessources res;

    load_patterns(res.regular_templates, res.partial_templates);
    init(&res, *mesh.get_grid_pointer());

    two_refinement_algorithm(*this, cellIdentifier, nb_levels);
  }
};





}


#endif