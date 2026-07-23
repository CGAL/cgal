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
#ifndef HEXMESHING_GET_CELLS_TO_REFINE_H
#define HEXMESHING_GET_CELLS_TO_REFINE_H

#include <CGAL/hexmeshing/Hexmeshing_two_refinement_mark_utils.h>
#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <iostream>

namespace CGAL::internal::Hexmeshing {
  /**
   * @brief Explores a face in a plane to identify marked nodes and gather adjacent faces
   *
   * This function performs a breadth-first exploration of faces in a plane, starting from
   * a given face. It identifies marked nodes (0-cells), gathers adjacent faces (2-cells),
   * and collects incident volumes (3-cells) that need refinement. The function:
   *   - Marks the current face as explored to avoid revisiting
   *   - Initializes the face's template_id to 0
   *   - Iterates through all edges of the face
   *   - Marks identified nodes that haven't been marked yet
   *   - Increments the face's template_id for each identified node
   *   - Explores unvisited edges and finds adjacent faces on the same plane
   *   - Adds adjacent faces to the exploration queue
   *   - Collects additional volumes that are incident to the explored edges
   *   - Adds the face to the plane's face list if it has any identified nodes
   *
   * This function is part of the plane exploration process that identifies which faces
   * and volumes need to be refined based on the identified nodes.
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data to store marked nodes, faces, and additional volumes
   * @param queue Queue of faces to explore (breadth-first traversal)
   * @param face The face to explore
   * @param explored_edge Mark to track explored edges
   * @param explored_face Mark to track explored faces
   */
  template <typename HexData>
  void explore_face_of_plane(HexData& hdata, RefinementData& rdata, std::queue<Dart_handle>& queue,
        Dart_handle face, size_type explored_edge, size_type explored_face) {
    LCC& lcc = hdata.lcc;

    auto& face_attr = lcc.attribute<2>(face)->info();

    if (!lcc.is_whole_cell_unmarked<2>(face, explored_face)) return ;
    lcc.mark_cell<2>(face, explored_face);

    face_attr.template_id = 0;

    // In multi threaded code, check if both volumes are owned, otherwise skip
    bool is_markable = true;

    // Might not be needed,

    // if constexpr (std::is_same_v<HexData, ProcessData>){
    //   auto& front_vol = lcc.attribute<3>(face)->info();
    //   auto back_vol_attr = lcc.attribute<3>(lcc.beta<3>(face));
    //   is_markable = front_vol.owned or (back_vol_attr != nullptr && back_vol_attr->info().owned);
    // }

    auto edges = lcc.darts_of_cell<2,1>(face);
    assert(edges.size() == 4);
    // Add neighboring faces
    for (auto dit = edges.begin(), dend = edges.end(); dit != dend; dit++){
      bool edge_explored = lcc.is_whole_cell_marked<1>(dit, explored_edge);
      bool marked = lcc.is_marked(dit, hdata.template_mark);
      bool identified = lcc.is_marked(dit, hdata.identified_mark);

      if (is_markable && !marked && identified){
        lcc.mark_cell<0>(dit, hdata.template_mark);
        rdata.marked_nodes.push_back(dit);
      }

      if (identified){
        face_attr.template_id++;
      }

      if (!edge_explored) {
        lcc.mark_cell<1>(dit, explored_edge);

        // Also add incident volumes while iterating
        Dart_handle other_face = __adjacent_face_on_plane(lcc, rdata.iteration, dit, rdata.additionnal_volumes_found);
        // Also do the other way around to catch volumes in the opposite dir
        if (!lcc.is_free<3>(dit)){
          Dart_handle other_face2 = __adjacent_face_on_plane(lcc, rdata.iteration, lcc.beta(dit, 3), rdata.additionnal_volumes_found);
          if (other_face == lcc.null_dart_descriptor) other_face  = other_face2;
        }

        if (other_face != lcc.null_dart_descriptor)
          queue.push(other_face);
      }
    }


    if (face_attr.template_id > 0) {
      rdata.faces_of_plane.push_back(face);
    }
  }

  /**
   * @brief Marks faces for propagation based on the template pattern of the given face
   *
   * This function identifies which faces need to be marked for propagation based on the
   * template_id of the given face. It handles two cases:
   * 1. For template_id=1: Finds the single marked node and marks three adjacent faces
   * 2. For template_id=2: Finds the edge with both endpoints marked and marks two adjacent faces
   *
   * The function locates the marked edge(s) on the face and then marks the appropriate
   * half-faces for propagation using the propagation_face_mark. This ensures that
   * refinement patterns can properly propagate to adjacent volumes.
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param face A dart handle representing the face to analyze for propagation marking
   * @param face_attr The face attribute containing the template_id information
   */
  template <typename HexData>
  void mark_template_for_propagation(HexData &hdata, Dart_handle face, DartInfo::FaceAttrValue &face_attr)
  {
    LCC& lcc = hdata.lcc;
    Dart_handle marked_edge;
    auto nodes = lcc.darts_of_cell<2, 0>(face);

    bool found = false;

    if (face_attr.template_id == 1)
    {
      for (auto it = nodes.begin(); it != nodes.end() && !found; it++)
      {
        if (lcc.is_marked(it, hdata.template_mark))
        {
          found = true;
          marked_edge = it;
        }
      }
    }

    else if (face_attr.template_id == 2)
    {
      for (auto it = nodes.begin(); it != nodes.end() && !found; it++)
      {
        if (lcc.is_marked(it, hdata.template_mark) && lcc.is_marked(lcc.other_extremity(it), hdata.template_mark))
        {
          found = true;
          marked_edge = it;
        }
      }
    }

    assert(found);

    Dart_handle face1 = lcc.beta(marked_edge, 2, 1, 1, 2);
    Dart_handle face2 = lcc.beta(face1, 0, 2);
    Dart_handle face3 = lcc.beta(face1, 0, 0, 2);

    mark_half_face_unchecked(lcc, face1, hdata.propagation_face_mark);
    mark_half_face_unchecked(lcc, face3, hdata.propagation_face_mark);

    if (face_attr.template_id == 1)
    {
      mark_half_face_unchecked(lcc, face2, hdata.propagation_face_mark);
    }
  }

  /**
   * @brief Propagates refinement from a marked face to adjacent volumes and faces
   *
   * This function handles the propagation of refinement patterns from a marked face
   * to its adjacent volumes and faces. It performs the following operations:
   * 1. Marks the current volume for refinement and sets its iteration
   * 2. Identifies the back face and its associated volume for refinement
   * 3. Adds the back volume to the refinement list if not already processed
   * 4. Marks the back face for refinement if not already explored
   * 5. Iterates through all edges of the back volume face
   * 6. Marks unmarked nodes and adds incident faces to the refinement list
   *
   * This function ensures that refinement patterns properly propagate through the
   * mesh structure, maintaining consistency across adjacent volumes and faces.
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data containing the current iteration and collections
   * @param face A dart handle representing the face to propagate from
   * @param explored_face_mark Mark used to track already explored faces
   */
  template <typename HexData>
  void propagate_face(HexData &hdata, RefinementData &rdata, const Dart_handle &face, size_type explored_face_mark)
  {
    LCC& lcc = hdata.lcc;


    assert(lcc.attribute<3>(face) != nullptr);
    auto &vol_attr = get_or_create_refinement_volume(lcc, face)->info();
    vol_attr.iteration = rdata.iteration;

    // /!\ ALSO: Mark add the hex attached to the back_face for refinement.
    // It is garanted that faces and volumes added will be unique in our array
    // Because that back_hex is only accessible within the template itself, and not
    // accessible from the plane.
    Dart_handle back_face = lcc.beta(face, 2, 1, 1, 2);
    Dart_handle back_volume_face = lcc.beta(back_face, 3);
    assert(back_volume_face != nullptr);
    auto &back_vol_attr = get_or_create_refinement_volume(lcc, back_volume_face)->info();

    if (back_vol_attr.iteration != rdata.iteration){
      rdata.volumes_to_refine.push_back(back_volume_face);
      back_vol_attr.iteration = rdata.iteration;
    }

    // Also add the neighboring faces to that 4-template face forrefinement
    if (!lcc.is_marked(back_volume_face, explored_face_mark)){
      rdata.faces_to_refine.push_back(back_volume_face);
      mark_face_unchecked(lcc, back_volume_face, explored_face_mark);
    }

    auto edges = lcc.darts_of_cell<2, 1>(back_volume_face);

    for (auto it = edges.begin(); it != edges.end(); it++)
    {
      auto top_face = lcc.beta(it, 2);

      if (!lcc.is_marked(it, hdata.template_mark)){
        lcc.mark_cell<0>(it, hdata.template_mark);
        rdata.marked_nodes.push_back(it);
      }

      if (!lcc.is_whole_cell_marked<2>(top_face, explored_face_mark)){
        rdata.faces_to_refine.push_back(top_face);
        mark_face_unchecked(lcc, top_face, explored_face_mark);
      }
    }
  }

  /**
   * @brief Executes the propagation stage of the refinement algorithm
   *
   * This function handles the propagation of refinement patterns across the mesh.
   * It consists of two main phases:
   *
   * Phase 1 (for iterations > 0): Propagates from 4-template faces
   * - Iterates through all faces in the current plane
   * - For faces with template_id=4, checks if they are marked for propagation
   * - Calls propagate_face for both the face and its opposite face (if it exists)
   *
   * Phase 2 (for iterations < 2): Marks faces for future propagation
   * - Iterates through all faces in the current plane
   * - For faces with template_id=1 or 2, calls mark_template_for_propagation
   * - Handles both the face and its opposite face (if it exists)
   *
   * The function skips propagation entirely for iteration 0 and stops marking
   * for propagation after iteration 2, as the refinement patterns are fully
   * established by that point.
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data containing the current iteration and faces of the plane
   * @param explored_face_mark Mark used to track already explored faces
   */
  template <typename HexData>
  void propagation_stage(HexData &hdata, RefinementData &rdata, size_type explored_face_mark)
  {
    LCC& lcc = hdata.lcc;

    bool skip_propagation = rdata.iteration == 0;
    int propagated_count = 0;
    int marked_for_prop_count = 0;

    if (!skip_propagation) for (Dart_handle face : rdata.faces_of_plane){
      auto& face_attr = lcc.attribute<2>(face)->info();

      if (face_attr.template_id == 4){
        if (is_half_face_marked(lcc, face, hdata.propagation_face_mark)){
          propagate_face(hdata, rdata, face, explored_face_mark);
          propagated_count++;
        }

        Dart_handle other_face = lcc.beta(face, 3);
        if (!lcc.is_free<3>(face) && is_half_face_marked(lcc, other_face, hdata.propagation_face_mark)) {
          propagate_face(hdata, rdata, other_face, explored_face_mark);
          propagated_count++;
        }
      }
    }

    if (rdata.iteration >= 2) return;

    for (Dart_handle face : rdata.faces_of_plane){
      auto& face_attr = lcc.attribute<2>(face)->info();
      if (face_attr.template_id < 1 or face_attr.template_id > 2) continue;

      mark_template_for_propagation(hdata, face, face_attr);
      marked_for_prop_count++;

      if (!lcc.is_free<3>(face)){
        mark_template_for_propagation(hdata, lcc.beta(face, 3), face_attr);
        marked_for_prop_count++;
      }
    }

    std::cout << "Number of faces propagated : " << propagated_count << std::endl;
    std::cout << "Number of faces marked for propagation : " << marked_for_prop_count << std::endl;
  }

  /**
   * @brief Collects cells that need refinement from the current plane
   *
   * This function processes all faces in the current refinement plane and identifies
   * which cells (faces and volumes) need to be refined. It performs the following operations:
   *
   * 1. For each face in the plane, iterates through all its edges
   * 2. For edges with marked nodes, finds incident faces normal to the plane
   * 3. Adds unmarked incident faces to the faces_to_refine collection
   * 4. Processes the volumes associated with each face:
   *    - Creates or updates volume attributes with the current iteration
   *    - Adds volumes to volumes_to_refine or partial_templates_to_refine based on template_id
   *    - Handles both the current face's volume and its opposite volume (if it exists)
   *
   * The function ensures that all cells affected by the refinement patterns are properly
   * identified and categorized for subsequent refinement operations.
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data containing the faces of the plane and collections to populate
   * @param explored_faces Mark used to track already explored faces to avoid duplicates
   */
  template <typename HexData>
  void get_cells_to_refine_from_plane(HexData &hdata, RefinementData &rdata, size_type explored_faces)
  {
    LCC& lcc = hdata.lcc;

    // Iterate over all faces of even planes
    for (Dart_handle face : rdata.faces_of_plane)
    {
      auto edges = lcc.darts_of_cell<2, 0>(face);

      for (auto edge = edges.begin(); edge != edges.end(); edge++)
      {
        // get incident faces to marked node to be refined
        // Incident faces normal to the plane
        // We also need to prevent adding twice the faces by marking them

        if (!lcc.is_marked(edge, hdata.template_mark) && !lcc.is_marked(lcc.other_extremity(edge), hdata.template_mark))
          continue;

        auto top_face_1 = lcc.beta(edge, 2);
        auto top_face_2 = lcc.beta(edge, 3, 2);

        if (!lcc.is_whole_cell_marked<2>(top_face_1, explored_faces))
        {
          rdata.faces_to_refine.push_back(top_face_1);
          mark_face_unchecked(lcc, top_face_1, explored_faces);
        }

        if (top_face_2 != lcc.null_dart_descriptor && !lcc.is_whole_cell_marked<2>(top_face_2, explored_faces))
        {
          rdata.faces_to_refine.push_back(top_face_2);
          mark_face_unchecked(lcc, top_face_2, explored_faces);
        }
      }
      // Also add the adjacent volumes if there is atleast one marked node
      // Because we are on odd/even layers, we can't accidently add twice a volume

      auto &face_attr = lcc.attribute<2>(face)->info();
      auto &vol_attr = get_or_create_refinement_volume(lcc, face)->info();

      bool vol_processed = vol_attr.iteration == rdata.iteration;

      // Two faces point to the same volume. We are on a merged cell (caused by 3 templates)
      // if (vol_processed){
      //   mark_all_0_cells<3>(lcc, face, hdata.template_mark);
      // }

      // Create first volume attr and push back the volume
      if (!vol_processed) {
        vol_attr.iteration = rdata.iteration;

        if (face_attr.template_id != 3){
          rdata.volumes_to_refine.push_back(face);
        }
        else
          rdata.partial_templates_to_refine.push_back(face); // Both volumes are treated if it is 3 template
      }

      // Create second volume attr and push the face from the other volume
      if (!lcc.is_free<3>(face))
      {
        Dart_handle other_face = lcc.beta(face, 3);
        auto &other_vol_attr = get_or_create_refinement_volume(lcc, other_face)->info();
        bool other_vol_processed = other_vol_attr.iteration == rdata.iteration;

        // if (other_vol_processed)
          // mark_all_0_cells<3>(lcc, other_face, hdata.template_mark);
        if (!other_vol_processed) {
          other_vol_attr.iteration = rdata.iteration;

          if (face_attr.template_id != 3){
            rdata.volumes_to_refine.push_back(other_face);
          }
        }
      }
    }
  }

  /**
   * @brief Collects cells that need refinement from additional volumes discovered during plane exploration
   *
   * This function processes additional volumes that were discovered during the plane exploration
   * phase but are not directly accessible from the current plane. These volumes typically
   * represent cells that are inside the refinement domain but not directly connected to
   * the current refinement plane.
   *
   * For each additional volume, the function:
   * 1. Checks if either endpoint of the volume's edge is marked
   * 2. If marked nodes are found, adds the volume to the refinement list
   * 3. Identifies the four adjacent faces of the volume
   * 4. Adds unmarked adjacent faces to the faces_to_refine collection
   * 5. Handles face marking based on which specific nodes are marked
   *
   * The function ensures that all volumes and faces affected by the refinement patterns
   * are properly identified, even if they are not directly accessible from the main plane.
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and marks
   * @param rdata Refinement data containing the additional volumes and collections to populate
   * @param explored_face Mark used to track already explored faces to avoid duplicates
   */
  template <typename HexData>
  void get_cells_to_refine_from_additionnal_volumes(HexData &hdata, RefinementData &rdata, size_type explored_face)
  {
    LCC& lcc = hdata.lcc;

    // No additionnal volumes should be found on the first iteration
    assert(rdata.iteration != 0 || rdata.iteration == 0 && rdata.additionnal_volumes_found.size() == 0);

    for (Dart_handle initial_edge : rdata.additionnal_volumes_found)
    {
      auto &vol_attr = lcc.attribute<3>(initial_edge)->info();

      bool node_1_marked = lcc.is_marked(initial_edge, hdata.template_mark);
      bool node_2_marked = lcc.is_marked(lcc.other_extremity(initial_edge), hdata.template_mark);

      Dart_handle adjacent_faces[] = {
        initial_edge,
        lcc.beta(initial_edge, 2),

        lcc.beta(initial_edge, 0, 2),
        lcc.beta(initial_edge, 1, 2)
      };

      if (node_1_marked || node_2_marked)
      {
        if (vol_attr.iteration != rdata.iteration)
          rdata.volumes_to_refine.push_back(initial_edge);

        vol_attr.iteration = rdata.iteration;

        if (!lcc.is_whole_cell_marked<2>(adjacent_faces[0], explored_face))
        {
          mark_face_unchecked(lcc, adjacent_faces[0], explored_face);
          rdata.faces_to_refine.push_back(adjacent_faces[0]);
        }

        if (!lcc.is_whole_cell_marked<2>(adjacent_faces[1], explored_face))
        {
          mark_face_unchecked(lcc, adjacent_faces[1], explored_face);
          rdata.faces_to_refine.push_back(adjacent_faces[1]);
        }
      }

      if (node_1_marked && !lcc.is_whole_cell_marked<2>(adjacent_faces[2], explored_face))
      {
        mark_face_unchecked(lcc, adjacent_faces[2], explored_face);
        rdata.faces_to_refine.push_back(adjacent_faces[2]);
      }

      if (node_2_marked && !lcc.is_whole_cell_marked<2>(adjacent_faces[3], explored_face))
      {
        mark_face_unchecked(lcc, adjacent_faces[3], explored_face);
        rdata.faces_to_refine.push_back(adjacent_faces[3]);
      }
    }
  }

  /**
   * @brief Main function to identify and collect all cells that need refinement for a given plane direction
   *
   * This function orchestrates the complete process of identifying cells that need refinement
   * for a specific plane direction (X, Y, or Z). It performs the following sequence of operations:
   *
   * 1. **Plane Exploration**: Explores all even planes in the specified direction using breadth-first search
   * 2. **Validation**: Validates that all faces in the plane have consistent template IDs
   * 3. **Conflict Resolution**: Fixes impossible template cases (diagonal 2-templates and adjacent 3-templates)
   * 4. **Node Communication**: Communicates marked nodes between threads (in parallel version)
   * 5. **Propagation**: Executes the propagation stage to mark faces for future propagation
   * 6. **Cell Collection**: Collects cells to refine from both the plane and additional volumes
   *
   * The function manages marks for tracking explored edges and faces, and ensures proper cleanup
   * of these marks at the end. This is the central coordination function for the refinement process.
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and plane information
   * @param rdata Refinement data to be populated with cells that need refinement
   * @param iterationPlane The plane direction (X, Y, or Z) for the current refinement iteration
   */
  template <typename HexData>
  void get_cells_to_refine(HexData &hdata, RefinementData &rdata, PlaneNormal iterationPlane)
  {
    using PlaneCC = std::vector<Dart_handle>;
    using PlaneSet = std::vector<PlaneCC>;

    LCC& lcc = hdata.lcc;
    rdata.iteration = iterationPlane;

    size_type explored_edge = lcc.get_new_mark();
    size_type explored_face = lcc.get_new_mark();

    PlaneSet& plane_set = hdata.first_face_of_planes[iterationPlane];

    // Explore all even planes
    for (int i = 1; i < plane_set.size(); i += 2) {
      std::queue<Dart_handle> to_explore;

      for (auto start : plane_set[i])
        to_explore.push(start);

      while (!to_explore.empty()) {
        Dart_handle front = to_explore.front();
        to_explore.pop();
        explore_face_of_plane(hdata, rdata, to_explore, front, explored_edge, explored_face);
      }
    }

    assert_faces_of_plane_valid(hdata, rdata);

    // Communicate found marked nodes,
    // fix procedure and propagation procedure will produce on two threads the same result

    fix_impossible_cases(hdata, rdata);
    assert_faces_of_plane_valid(hdata, rdata);

    thread_communicate_marked_nodes(hdata, rdata, explored_edge);
    propagation_stage(hdata, rdata, explored_face);

    get_cells_to_refine_from_plane(hdata, rdata,  explored_face);
    get_cells_to_refine_from_additionnal_volumes(hdata, rdata, explored_face);

    lcc.free_mark(explored_edge);
    lcc.free_mark(explored_face);

  }
}





#endif