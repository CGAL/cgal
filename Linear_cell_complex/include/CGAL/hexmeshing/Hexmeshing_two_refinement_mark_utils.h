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
#ifndef CGAL_HEXMESHING_TWO_REFINEMENT_MARK_UTILS_H
#define CGAL_HEXMESHING_TWO_REFINEMENT_MARK_UTILS_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_refinement_data.h>
#include <CGAL/hexmeshing/Hexmeshing_two_refinement_utils.h>
#include <cassert>

namespace CGAL::internal::Hexmeshing {
  /**
   * @brief Marks all k-cells that are part of a given i-cell
   *
   * This　function iterates over all darts of an i-dimensional cell and marks
   * all k-dimensional cells that are part of it.
   *
   * @tparam i The dimension of the source cell
   * @tparam k The dimension of the cells to mark
   * @param lcc The Linear Cell Complex
   * @param dart A dart handle representing the i-cell
   * @param mark The mark to apply to the k-cells
   */
  template <uint i, uint k>
  void mark_k_cells_of_i_cell(LCC& lcc, Dart_handle dart, size_type mark){
    auto iterator = lcc.darts_of_cell<i, 0>(dart);
    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++){
      if (!lcc.is_marked(dit, mark)) {
        lcc.mark_cell<k>(dit, mark);
      }
    }
  }

  /**
   * @brief Marks a face without checking its current mark status
   *
   * This function marks all darts of a 2-dimensional cell (face) with the specified mark.
   * Unlike safe marking functions, this function does not verify if the face is already marked.
   *
   * @param lcc The Linear Cell Complex
   * @param dart A dart handle representing the face to mark
   * @param mark The mark to apply to the face
   */
  void mark_face_unchecked(LCC& lcc, Dart_handle dart, size_type mark){
    auto iterator = lcc.darts_of_cell<2, 1>(dart);
    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++) {
      lcc.mark(dit, mark);

      if (!lcc.is_free<3>(dit))
        lcc.mark(lcc.beta<3>(dit), mark);
    }

    assert(lcc.is_whole_cell_marked<2>(dart, mark));
  }

  /**
   * @brief Checks if a half-face is marked with the specified mark
   *
   * This function checks if all darts of a half-face are marked with the given mark.
   * A half-face consists of the 4 darts that can be reached from the specified dart
   * without using beta<3> operations, out of the total 8 darts that include beta<3>.
   *
   * @param lcc The Linear Cell Complex
   * @param dart A dart handle representing the half-face to check
   * @param mark The mark to check for
   * @return true if all darts of the half-face are marked, false otherwise
   */
  bool is_half_face_marked(const LCC& lcc, LCC::Dart_const_handle dart, size_type mark ){
    auto iterator = lcc.darts_of_cell<2, 1>(dart);
    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++) {
      if (!lcc.is_marked(dit, mark)) return false;
    }

    return true;
  }

  /**
   * @brief Marks a half-face without checking its current mark status
   *
   * This function marks all darts of a half-face with the specified mark.
   * A half-face consists of the 4 darts that can be reached from the specified dart
   * without using beta<3> operations, out of the total 8 darts that include beta<3>.
   * Unlike safe marking functions, this function does not verify if the half-face is already marked.
   *
   * @param lcc The Linear Cell Complex
   * @param dart A dart handle representing the half-face to mark
   * @param mark The mark to apply to the half-face
   */
  void mark_half_face_unchecked(LCC& lcc, Dart_handle dart, size_type mark){
    auto iterator = lcc.darts_of_cell<2, 1>(dart);

    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++) {
      lcc.mark(dit, mark);
    }

    assert(is_half_face_marked(lcc, dart, mark));
  }

  /**
   * @brief Expands the marking of a 0-cell (node) to affect surrounding faces
   *
   * This function handles the propagation of node marking by finding all faces
   * around the marked node and updating their template IDs. When a node is marked,
   * all faces that share this node need to have their template_id incremented.
   * Faces that transition from template_id 0 to 1 are added to the refinement list,
   * and faces with template_id 2 or 3 are added to the faces_to_check queue for
   * further processing.
   *
   * @param lcc The Linear Cell Complex
   * @param rdata Refinement data containing the current iteration direction and collections
   * @param faces_to_check Queue of faces that need to be checked for template conflicts
   * @param edge A dart handle representing the marked node (0-cell)
   */
  void __expand_0_cell_marking(LCC &lcc, RefinementData &rdata, std::queue<Dart_handle>& faces_to_check, Dart_handle &edge) {
    auto faces = plane_faces_around_node(lcc, rdata, edge);
    int s = faces.size();

    for (Dart_handle face : faces){
      assert( lcc.attribute<2>(face) != nullptr);
      auto& face_attr =  lcc.attribute<2>(face)->info();

      // If the face didn't have any template before, it will have one, so add it in faces to refine
      if (face_attr.template_id == 0) {
        rdata.faces_of_plane.push_back(face);
      }

      face_attr.template_id++;
      assert(face_attr.template_id <= 4);

      if (face_attr.template_id == 2 || face_attr.template_id == 3)
        faces_to_check.push(face);
    }
  }

  /**
   * @brief Fixes diagonal marking patterns on faces with template_id=2 by converting them to template_id=4
   *
   * This function handles the specific case where a face has template_id=2 with marks placed
   * diagonally (e.g., nodes 1 and 3 are marked but node 2 is not). In such cases, the
   * refinement pattern becomes impossible to resolve, so the function marks all unmarked
   * nodes on the face to create a complete 4-template pattern.
   *
   * The function examines the face's edges in sequence and detects diagonal marking patterns
   * where marked nodes are not consecutively chained. If such a pattern is found, it marks
   * all remaining unmarked nodes and sets the face's template_id to 4.
   *
   * @pre The given face must have template_id=2
   * @param lcc the linear cell complex
   * @param rdata Refinement data containing the current iteration direction
   * @param template_mark HexMeshingData template_mark
   * @param faces_to_check Queue of faces that need to be checked for template conflicts
   * @param face A dart handle representing the face with template_id=2 to check for diagonal marking
   * @return true if a fix was applied (diagonal pattern was converted to template_id=4), false if the marks were already consecutive
   */
  bool fix_mark_connectivity(LCC &lcc, RefinementData &rdata, size_type template_mark, std::queue<Dart_handle>& faces_to_check, Dart_handle face){
    auto edges = lcc.darts_of_cell<2, 1>(face);

    bool connected = true;
    Dart_handle edge1 = face;
    if (lcc.is_marked(edge1, template_mark)
      && !lcc.is_marked(lcc.beta(edge1, 1), template_mark)
      &&  lcc.is_marked(lcc.beta(edge1, 1, 1), template_mark))
    { connected = false;}

    Dart_handle edge2 = lcc.beta(face, 1);
    if (lcc.is_marked(edge2, template_mark)
      && !lcc.is_marked(lcc.beta(edge2, 1), template_mark)
      &&  lcc.is_marked(lcc.beta(edge2, 1, 1), template_mark))
    { connected = false;}

    if (connected) return false;

    // mark the face
    for (auto edge = edges.begin(); edge != edges.end(); edge++) {
      if (lcc.is_marked(edge, template_mark)) continue;

      lcc.mark_cell<0>(edge, template_mark);
      rdata.marked_nodes.push_back(edge);
      __expand_0_cell_marking(lcc, rdata, faces_to_check, edge);
    }

    auto &face_attr = lcc.attribute<2>(face)->info();
    face_attr.template_id = 4;

    return true;

  }

  /**
   * @brief Marks all nodes (0-cells) that belong to identified volumes for refinement
   *
   * This function iterates through all volume attributes in the Linear Cell Complex
   * and marks all nodes (vertices) that belong to volumes that have been identified
   * for refinement. The marking is done using the identified_mark to track which
   * nodes should be considered during the refinement process.
   *
   * The function performs the following operations:
   *
   * 1. **Volume Attribute Iteration**: Iterates through all volume attributes (3-cells)
   *    in the Linear Cell Complex using `lcc.attributes<3>()`
   *
   * 2. **Identification Check**: For each volume, checks if:
   *    - The volume's type is greater than `VolumeType::NONE` (i.e., it has been
   *      identified for refinement, is in expansion, or is marked for refinement)
   *    - The volume is owned by the current process (`owned` flag is true)
   *
   * 3. **Node Marking**: For volumes that meet the criteria, calls `mark_k_cells_of_i_cell<3, 0>`
   *    to mark all nodes (0-cells) that belong to the volume with the `identified_mark`
   *
   * This function is called at the beginning of each plane iteration in the refinement
   * algorithm to ensure that all nodes belonging to identified volumes are properly
   * marked before the refinement process begins. This marking is essential for the
   * subsequent template identification and refinement operations.
   *
   * @param lcc the linear cell complex
   * @param identified_mark HexMeshingData identified_mark
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   */
  void mark_identified_cells_from_3_attrs(LCC& lcc, size_type identified_mark) {
    auto& attributes = lcc.attributes<3>();

    for (auto it = attributes.begin(), end = attributes.end(); it != end; it++){
      if (it->info().type > VolumeType::NONE && it->info().owned){
        mark_k_cells_of_i_cell<3, 0>(lcc, it->dart(), identified_mark);
      }
    }
  }
}

#endif // CGAL_HEXMESHING_TWO_REFINEMENT_MARK_UTILS_H