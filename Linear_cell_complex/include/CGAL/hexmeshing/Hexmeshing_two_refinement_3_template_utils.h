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
#ifndef CGAL_HEXMESHING_TWO_REFINEMENT_3_TEMPLATE_UTILS_H
#define CGAL_HEXMESHING_TWO_REFINEMENT_3_TEMPLATE_UTILS_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_two_refinement_mark_utils.h>
#include <CGAL/hexmeshing/Hexmeshing_threads_utils.h>
#include <cassert>

namespace CGAL::internal::Hexmeshing {
  /**
   * @brief Finds the origin dart of a 3-template pattern
   *
   * This function searches for the origin dart of a 3-template pattern by traversing
   * through the marked faces.
   *
   * @param lcc The Linear Cell Complex
   * @param marked_face A dart handle representing the marked face to start the search from
   * @param template_mark The mark used to identify template faces
   * @return A dart handle pointing to the origin of the 3-template pattern, or null if not found
   */
  Dart_handle find_3_template_origin(LCC& lcc, Dart_handle marked_face, size_type template_mark) {

    auto _edges_count = lcc.darts_of_cell<2,1>(marked_face).size();
    assert(_edges_count == 6);

    Dart_handle dart = marked_face;

    // Get the origin dart : Find the two unmarked node on the face
    // since the 3 template is created by adjacent two 2-templates
    bool found = false;
    for (int i = 0; i < 6; i++)
    {
      if (!lcc.is_marked(dart, template_mark)
        && !lcc.is_marked(lcc.beta(dart, 1), template_mark)
        && !lcc.is_marked(lcc.beta(dart, 1, 1), template_mark))
      {
        assert(lcc.is_marked(lcc.beta(dart, 1, 1, 1), template_mark));
        assert(lcc.is_marked(lcc.beta(dart, 1, 1, 1, 1), template_mark));
        assert(lcc.is_marked(lcc.beta(dart, 1, 1, 1, 1, 1), template_mark));
        assert(lcc.beta(dart, 1, 1, 1, 1, 1, 1) == dart);
        found = true;
        break;
      }

      dart = lcc.beta(dart, 1);
    }

    assert(found);

    return dart;
  }

  /**
   * @brief Fixes adjacent 3-templates by transforming them into 4-templates when they share an unmarked node
   *
   * This function handles the case where two adjacent faces both have template_id 3 and share
   * an unmarked node. In such cases, the refinement pattern becomes impossible to resolve,
   * so the function marks the shared unmarked node to force both faces to become 4-templates.
   *
   * The function finds the unmarked edge on the given face, checks if the adjacent face
   * on the same plane also has template_id 3, and if so, marks the unmarked node to resolve
   * the conflict.
   *
   * @param lcc the linear cell complex
   * @param rdata Refinement data containing the current iteration direction
   * @param template_mark HexMeshingData template_mark
   * @param faces_to_check Queue of faces that need to be checked for template conflicts
   * @param face A dart handle representing the face to check for adjacent 3-template conflicts
   * @return true if a fix was applied (node was marked), false otherwise
   */
  bool fix_adjacent_3_templates(LCC &lcc, RefinementData &rdata, size_type template_mark, std::queue<Dart_handle>& faces_to_check, Dart_handle &face)
  {
    // Transform nearby 3 templates into two 4 templates
    // if they both share the unmarked node

    auto &face_attr = lcc.attribute<2>(face)->info();

    auto edges = lcc.darts_of_cell<2, 1>(face);

    // find the unmarked edge;
    bool found = false;
    Dart_handle edges_unmarked[2];

    for (auto it = edges.begin(); it != edges.end() && !found; it++)
    {
      if (!lcc.is_marked(it, template_mark)){
        found = true;
        edges_unmarked[0] = it;
      }
    }

    assert(found);
    edges_unmarked[1] = lcc.beta(edges_unmarked[0], 0);

    bool is_impossible = false;
    for (int i = 0; i < 2; i++){
      Dart_handle edge = edges_unmarked[i];
      Dart_handle other_face = adjacent_face_on_plane(lcc, rdata.iteration, edge);

      if (other_face == lcc.null_dart_descriptor) continue;

      auto other_face_handle = lcc.attribute<2>(other_face);

      assert(other_face_handle != nullptr);

      if (other_face_handle->info().template_id == 3){
        is_impossible = true;
        break;
      }

    }

    if (is_impossible) {
      lcc.mark_cell<0>(edges_unmarked[0], template_mark);
      rdata.marked_nodes.push_back(edges_unmarked[0]);
      __expand_0_cell_marking(lcc, rdata, faces_to_check, edges_unmarked[0]);

      return true;
    }

    return false;
  }

  /**
   * @brief Cleans up a 3-template by removing temporary elements and merging face attributes
   *
   * This function performs the cleanup phase after a 3-template refinement operation.
   * It removes temporary elements created during the refinement process and merges
   * the face attributes of the two neighboring faces that were split during the
   * 3-template creation.
   *
   * The function performs the following operations:
   * 1. Preserves face attributes before they are potentially destroyed
   * 2. Inserts barycenters in the lower edge and contracts them
   * 3. Contracts the two remaining 2-dart faces that were created during refinement
   * 4. Removes the temporary volume created by the partial template substitution
   * 5. Sews the two neighboring volumes together if both faces are valid
   * 6. Merges the plane information from the two original faces into a single face attribute
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex
   * @param origin_dart The origin dart of the 3-template volume to be cleaned up
   * @param upper_edge The upper edge of the 3-template
   * @param lower_edge The lower edge of the 3-template
   * @param face1 Reference to the first neighboring face (may be modified)
   * @param face2 Reference to the second neighboring face (may be modified)
   */
  template <typename HexData>
  void clean_up_3_template(HexData &hdata, const Dart_handle &origin_dart, const Dart_handle &upper_edge, const Dart_handle lower_edge, Dart_handle &face1, Dart_handle &face2)
  {
    LCC& lcc = hdata.lcc;

    // Ne pas maintenir des refs sur les attributs, ils vont disparaitre
    DartInfo::FaceAttrValue face1_attr, face2_attr;
    if (lcc.attribute<2>(face1) != nullptr) face1_attr = lcc.attribute<2>(face1)->info();
    if (lcc.attribute<2>(face2) != nullptr) face2_attr = lcc.attribute<2>(face2)->info();

    Dart_handle lower_mid_1 = lcc.insert_barycenter_in_cell<1>(lower_edge);
    Dart_handle lower_mid_2 = lcc.beta(lower_mid_1, 2, 1);

    thread_join_3_template_vertex__pairpair(hdata, lower_edge);

    lcc.contract_cell<1>(lower_mid_1);
    lcc.contract_cell<1>(lower_mid_2);

    // Contract the two remaining 2-darts faces
    Dart_handle face_to_remove_1 = lcc.beta(upper_edge, 2, 1);
    Dart_handle face_to_remove_2 = lcc.beta(upper_edge, 2, 1, 1, 2, 1, 2);

    assert(lcc.darts_of_orbit<1>(face_to_remove_1).size() == 2);
    assert(lcc.darts_of_orbit<1>(face_to_remove_2).size() == 2);

    // Contract the two remaining 2-darts faces
    lcc.contract_cell<2>(face_to_remove_1);
    lcc.contract_cell<2>(face_to_remove_2);

    // Remove the created volume from the partial query replace and sew the two neighboring volumes
    lcc.remove_cell<3>(origin_dart);

    if (face1 != lcc.null_dart_descriptor && face2 != lcc.null_dart_descriptor){
      lcc.sew<3>(face1, face2);
      assert(lcc.attribute<2>(face1) == lcc.attribute<2>(face2));
    }

    // Requires at least one face for merging/adding 2-attributes
    std::bitset<3> merged_planes = face1_attr.plane | face2_attr.plane;
    auto merge_face = face1 != lcc.null_dart_descriptor ? face1 : face2;

    assert(merge_face != lcc.null_dart_descriptor);

    // Merge the two previous face attributes bitsets
    auto& merge_face_attr = get_or_create_attr<2>(lcc, merge_face)->info();
    // TODO don't entirely reset...
    merge_face_attr = DartInfo::FaceAttrValue();
    merge_face_attr.plane = merged_planes;
  }

  /**
   * @brief Refines all 3-template patterns in the current refinement data
   *
   * This function processes all partial templates (3-templates) collected in the
   * refinement data and applies the necessary topological and combinatorial operations
   * to convert each 3-template into two regular hexahedral volumes. For each 3-template:
   *   - The function identifies the origin dart and its adjacent volume
   *   - It locates the relevant edges and faces for both volumes
   *   - It inserts barycenters into the upper edge and contracts the resulting edges
   *   - It calls clean_up_3_template to remove temporary elements and merge face attributes
   *   - The process is repeated for both the original and adjacent volumes
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex
   * @param rdata Refinement data containing the list of partial templates to refine
   * @return The total number of volumic substitutions performed (should be twice the number of 3-templates)
   */
  template <typename HexData>
  int refine_3_template(HexData &hdata, RefinementData& rdata)
  {
    int nb_sub = 0;
    for (auto origin_dart : rdata.partial_templates_to_refine){
      LCC& lcc = hdata.lcc;

      int nb_edges = lcc.darts_of_cell<2,1>(origin_dart).size();
      assert(nb_edges == 3);

      Dart_handle vol2_origin_dart = lcc.beta(origin_dart, 3);

      Dart_handle upper_edge = lcc.beta(origin_dart, 0); //?
      Dart_handle vol2_upper_edge = lcc.beta(upper_edge, 3);

      // Assert the origin dart have different directions on the adjacent volume
      // So assert that their opposite are equal
      assert(lcc.beta(origin_dart, 1) == lcc.beta(vol2_origin_dart, 0, 3));

      // Face of the two neighboring volumes to the created volume
      Dart_handle face1 = lcc.beta(origin_dart, 2, 3);
      Dart_handle face2 = lcc.beta(origin_dart, 1, 2, 3);
      Dart_handle lower_edge = lcc.beta(upper_edge, 2, 1, 1);

      Dart_handle vol2_face1 = lcc.beta(vol2_origin_dart, 2, 3);
      Dart_handle vol2_face2 = lcc.beta(vol2_origin_dart, 0, 2, 3); // 0 because opposite directions
      Dart_handle vol2_lower_edge = lcc.beta(vol2_upper_edge, 2, 1, 1);

      // Contract upper and lower edge into its barycenter

      Dart_handle upper_mid_1 = lcc.insert_barycenter_in_cell<1>(upper_edge);
      Dart_handle upper_mid_2 = lcc.beta(upper_mid_1, 2, 1);

      thread_join_3_template_vertex__pair(hdata, upper_edge);

      // contract the shared edge between the two volumes
      lcc.contract_cell<1>(upper_mid_1);
      lcc.contract_cell<1>(upper_mid_2);

      clean_up_3_template(hdata, origin_dart, upper_edge, lower_edge, face1, face2);
      clean_up_3_template(hdata, vol2_origin_dart, vol2_upper_edge, vol2_lower_edge, vol2_face1, vol2_face2);
      nb_sub += 2;
    }

    return nb_sub;
  }

  /**
   * @brief Refines all partial 3-template patterns using partial template substitution
   *
   * This function processes all faces marked as partial templates (3-templates) in the
   * refinement data and applies partial template substitution to each. For each such face:
   *   - The origin dart of the 3-template is found and used for substitution
   *   - The missing edge is created to enable the query_replace operation
   *   - The partial template substituter is used to replace the volume(s) with the correct pattern
   *   - Both the original and the adjacent 3-template volumes are processed
   *
   * This operation is necessary to handle cases where a regular template cannot be applied
   * due to the presence of a 3-template configuration.
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and external resources
   * @param rdata Refinement data containing the list of partial templates to refine
   */
  template <typename HexData>
  void refine_partial_templates(HexData& hdata, RefinementData& rdata)
  {
    LCC& lcc = hdata.lcc;

    // Partially refine the 3 templates
    for (auto& marked_face : rdata.partial_templates_to_refine){
      assert(lcc.attribute<2>(marked_face) != nullptr);
      assert(lcc.attribute<2>(marked_face)->info().template_id == 3);

      // Query replace with the partial 3-template, making it into two volumes
      Dart_handle origin_dart = find_3_template_origin(lcc, marked_face, hdata.template_mark);
      marked_face = origin_dart;

      Dart_handle upper_d1 = origin_dart;
      Dart_handle upper_d2 = lcc.beta(origin_dart, 1, 1);

      // Create the missing edge to perform the query_replace
      Dart_handle upper_edge = lcc.insert_cell_1_in_cell_2(upper_d1, upper_d2);

      size_type p = hdata.ext->partial_templates.query_replace_one_volume(lcc, marked_face, hdata.template_mark);
      assert(p == 0);

      // Also replace the other connected volume that is 3 template
      p = hdata.ext->partial_templates.query_replace_one_volume(lcc, lcc.beta(marked_face, 3), hdata.template_mark);
      assert(p == 0);
    }
  }

}


#endif // CGAL_HEXMESHING_TWO_REFINEMENT_3_TEMPLATE_UTILS_H