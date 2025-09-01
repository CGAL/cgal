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
#ifndef HEXMESHING_UTILS_H
#define HEXMESHING_UTILS_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_refinement_data.h>
#include <CGAL/hexmeshing/Hexmeshing_adjacent_face_on_plane.h>
#include <CGAL/hexmeshing/Hexmeshing_function_alias.h>
#include <CGAL/Union_find.h>
#include <queue>
#include <array>

#include <boost/container/static_vector.hpp>


namespace CGAL::internal::Hexmeshing {
  template <typename T>
  // TODO make union_find arg const
  std::vector<typename CGAL::Union_find<T>::handle> get_partitions(CGAL::Union_find<T>& union_find){
    using Handle = typename CGAL::Union_find<T>::handle;
    std::vector<Handle> result;

    if (union_find.number_of_sets() == 0) return result;

    for (auto it = union_find.begin(), end = union_find.end(); it != end; it++ ){
      if (it.ptr()->up == nullptr) {
        result.push_back(it);

        if (result.size() == union_find.number_of_sets())
          return result;
      }
    }

    CGAL_assertion_msg(false, "get_partions(Union_find<T>) function did not find all sets. Bug?");
    return result;
  };

  /**
   * @brief Internal function to iterate over faces and their edges in a plane
   *
   * This function performs a breadth-first traversal of faces in a plane, ensuring each face
   * and edge is processed exactly once. It uses marks to track visited elements and applies
   * the provided operations to faces and edges.
   *
   * @tparam FaceOp Type of the face operation functor
   * @tparam EdgeOp Type of the edge operation functor
   * @param lcc The Linear Cell Complex
   * @param to_explore Queue of dart handles to explore
   * @param face_operation Operation to apply to each face
   * @param edge_operation Operation to apply to each edge
   */
  template <typename FaceOp, typename EdgeOp>
  inline void __plane_for_each_face(LCC& lcc, std::queue<Dart_handle>& to_explore,
                                const FaceOp&& face_operation,
                                const EdgeOp&& edge_operation)
  {
    size_type face_mark = lcc.get_new_mark();
    size_type edge_mark = lcc.get_new_mark();

    while (!to_explore.empty()){
      Dart_handle face = to_explore.front();
      to_explore.pop();

      if (lcc.is_whole_cell_marked<2>(face, face_mark)) continue;
      lcc.mark_cell<2>(face, face_mark);

      auto edges = lcc.darts_of_cell<2,1>(face);
      face_operation(face, edges);

      for (auto it = edges.begin(), end = edges.end(); it != end; it++){
        if (lcc.is_whole_cell_marked<1>(it, edge_mark)) continue;
        lcc.mark_cell<1>(it, edge_mark);
        Dart_handle edge = edge_operation(it);
        if (edge != lcc.null_dart_descriptor)
          to_explore.push(edge);
      }
    }

    lcc.free_mark(face_mark);
    lcc.free_mark(edge_mark);
  }

  /**
   * @brief Retrieves all faces of a hexahedral cell
   *
   * This function returns an array containing dart handles for all six faces of a hexahedron.
   * The faces are returned in a specific order:
   * 1. The input face (vol)
   * 2-5. The faces adjacent to the edges of the input face
   * 6. The opposite face
   *
   * @param lcc The Linear Cell Complex
   * @param vol A dart handle representing the initial face of the hexahedron
   * @return std::array<Dart_handle, 6> Array of dart handles for all faces
   */
  std::array<Dart_handle, 6> faces_of_hex(LCC& lcc, Dart_handle vol){
    std::array<Dart_handle, 6> arr;
    int i = 0;

    arr[i++] = vol;

    auto edges = lcc.darts_of_cell<2,1>(vol);
    for (auto it = edges.begin(), end = edges.end(); it != end; it++){
        arr[i++] = lcc.beta(it, 2);
    }

    arr[i++] = lcc.beta(vol, 2, 1, 1, 2);

    return arr;
  }

  /**
   * @brief Iterates over all faces in a plane and applies specified operations
   *
   * This function performs a breadth-first traversal of faces in a plane, starting from
   * a given face. For each face and its edges, it applies the provided operation functors.
   * The function ensures that each face and edge is processed exactly once by using marks.
   *
   * @tparam FaceOp Type of the face operation functor
   * @tparam EdgeOp Type of the edge operation functor
   * @param lcc The Linear Cell Complex
   * @param start The starting face's dart handle
   * @param face_operation Operation to apply to each face
   * @param edge_operation Operation to apply to each edge
   */
  template <typename FaceOp, typename EdgeOp>
  void plane_for_each_face(LCC& lcc, std::vector<Dart_handle> starts,
                            const FaceOp&& face_operation,
                            const EdgeOp&& edge_operation)
  {
    std::queue<Dart_handle> to_explore;
    for (Dart_handle start : starts){
      to_explore.push(start);
    }

    __plane_for_each_face(lcc, to_explore,
      std::forward<const FaceOp>(face_operation),
      std::forward<const EdgeOp>(edge_operation));
  }

  /**
   * @brief Gets or creates an attribute for a cell in the Linear Cell Complex
   *
   * This template function retrieves an existing attribute for a cell, or creates a new one
   * if it doesn't exist. The attribute type and cell dimension are specified through template
   * parameters.
   *
   * @tparam i The dimension of the cell (0 for vertex, 1 for edge, 2 for face, 3 for volume)
   * @param lcc The Linear Cell Complex
   * @param dart A dart handle representing the cell
   * @return The attribute descriptor for the cell
   */
  template <uint i>
  typename LCC::Attribute_descriptor<i>::type get_or_create_attr(LCC& lcc, Dart_handle dart){
    auto attr = lcc.attribute<i>(dart);

    if (attr == nullptr){
      attr = lcc.create_attribute<i>();
      lcc.set_attribute<i>(dart, attr);
    }
    return attr;
  }

  /**
   * @brief Gets or creates a volume attribute for refinement operations
   *
   * This function retrieves an existing volume attribute or creates a new one if it doesn't exist.
   * The created volume attribute is initialized with default values suitable for refinement:
   * - type is set to VolumeType::REFINEMENT
   * - iteration is set to -1
   * - cc_id is set to the maximum possible value
   *
   * @param lcc The Linear Cell Complex
   * @param dart A dart handle representing the volume
   * @return The volume attribute descriptor
   */
  LCC::Attribute_descriptor<3>::type get_or_create_refinement_volume(LCC& lcc, Dart_handle dart){
    auto attr = get_or_create_attr<3>(lcc, dart);

    // Previously NONE volumes are tagged as refined
    if (attr->info().type == VolumeType::NONE)
      attr->info().type = VolumeType::REFINEMENT;

    return attr;

  }

  /**
   * @brief Retrieves all 26-connected neighboring cells of a given cell
   *
   * This function collects all cells that share a vertex, edge, or face with the given cell.
   * In a 3D grid, a cell can have up to 26 neighbors.
   *
   * @param lcc The Linear Cell Complex
   * @param dart A dart handle representing the central cell
   * @param include_self_vol Whether to include the input cell in the result (default: false)
   * @return A static vector containing dart handles for all neighboring cells (max size: 27)
   */
  boost::container::static_vector<Dart_handle, 27> cells_26_connectivity(LCC& lcc, Dart_handle dart, bool include_self_vol = false) {
    boost::container::static_vector<Dart_handle, 27> array;
    Dart_handle layers[3] = {
      dart,
      lcc.beta(dart, 3),
      lcc.beta(dart, 2, 1, 1, 2, 3)
    };

    for (int i = 0; i < 3; i++){
      Dart_handle mid_face_dart = layers[i];

      if (mid_face_dart == lcc.null_dart_descriptor) continue;

      if (include_self_vol || i != 0) array.push_back(mid_face_dart);

      auto edges = lcc.darts_of_cell<2, 1>(mid_face_dart);
      for (auto edge_it = edges.begin(), end = edges.end(); edge_it != end; edge_it++){

        if (!lcc.is_free<3>(lcc.beta(edge_it, 2))){
          Dart_handle other_face = lcc.beta(edge_it, 2, 3, 2);
          if (other_face != lcc.null_dart_descriptor)
            array.push_back(other_face);

          // TODO If nesting ..
          if (!lcc.is_free<3>(lcc.beta(other_face, 1, 2))){
            Dart_handle side_face = lcc.beta(other_face, 1, 2, 3, 2);
            if (side_face != lcc.null_dart_descriptor)
              array.push_back(side_face);
          }
        }

      }
    }

    return array;
  }

  /**
   * @brief Finds all faces around a node that lie on a specific plane
   *
   * This function collects all faces that are adjacent to a given node and lie on the specified plane.
   * There may be up to 8 faces adjacent to a vertex, where each dart handle belongs to a unique edge
   * and a unique face on the plane around the selected node. The function performs a traversal around
   * the node to gather all connected faces on the plane.
   *
   * @param lcc The Linear Cell Complex
   * @param rdata Refinement data containing the current iteration direction
   * @param node A dart handle representing the node to find faces around
   * @return A static vector containing dart handles for all faces around the node on the plane (max size: 8)
   */
  boost::container::static_vector<Dart_handle, 8>
  plane_faces_around_node(LCC& lcc,
                                    RefinementData &rdata,
                                    Dart_handle node){
    boost::container::static_vector<Dart_handle, 8> arr;

    // Add initial face
    arr.push_back(node);
    auto this_face_attr = lcc.attribute<2>(node);

    auto turn_around_edge = [&](Dart_handle edge){
      Dart_handle adjacent_face = adjacent_face_on_plane(lcc, rdata.iteration, edge);
      auto other_face_attr = lcc.attribute<2>(adjacent_face);

      while (adjacent_face != lcc.null_dart_descriptor && this_face_attr != other_face_attr){
        arr.push_back(adjacent_face);

        int f = lcc.belong_to_same_cell<0>(adjacent_face, node) ? 0 : 1;
        adjacent_face = adjacent_face_on_plane(lcc, rdata.iteration, lcc.beta(adjacent_face, f));
        if (adjacent_face != lcc.null_dart_descriptor){
          // Check if the dart is expectedly turning around the node
          assert(lcc.belong_to_same_cell<0>(adjacent_face, node) or lcc.belong_to_same_cell<0>(lcc.other_extremity(adjacent_face), node));
        }
        other_face_attr = lcc.attribute<2>(adjacent_face);
      } ;

      return adjacent_face;
    };

    // Turn around the first edge, gathering faces encountered
    Dart_handle end = turn_around_edge(node);

    // If the iteration doesn't loop back to the start, do it again on the opposite edge to that node
    if (lcc.attribute<2>(end) != this_face_attr){
      end = turn_around_edge(lcc.beta(node, 0));
    }

    // Check if all found faces are unique, if they are not, this means that the LCC was not properly refined
    auto assertion = [&](){
      std::unordered_set<DartInfo::FaceAttrValue*> faces_set;
      for (Dart_handle face : arr){
        auto f = lcc.attribute<2>(face);
        CGAL_postcondition_msg(f != nullptr, "plane_faces_around_node: returned array contains nullptr, Is the refinement correctly done ?");
        CGAL_postcondition_msg(true, "plane_faces_around_node: array contains a duplicate face, Is the refinement correctly done ?");
        assert(faces_set.count(&f->info()) == 0);
        faces_set.insert(&f->info());
      }
    };

    CGAL_postcondition_code(assertion());

    return arr;
  }

  /**
   * @brief Removes volumes from the Linear Cell Complex based on a trimming function
   *
   * This function iterates through all volumes in the Linear Cell Complex and removes
   * those that do not satisfy the criteria specified by the trimming function.
   *
   * The function performs the following operations:
   *
   * 1. **Volume Iteration**: Iterates through all volumes (3-cells) in the Linear Cell Complex
   *    using `lcc.one_dart_per_cell<3>()`
   *
   * 2. **Volume Evaluation**: For each volume, calls the provided trimming function
   *    `func(lcc, it)` to determine if the volume should be kept
   *
   * 3. **Volume Removal**: If the trimming function returns false (indicating the volume
   *    should be removed), calls `lcc.remove_cell<3>(it)` to remove the volume from the mesh
   *
   * This function is typically used to remove volumes that are outside the desired domain
   * or do not meet certain geometric or topological criteria. It is commonly called after
   * the refinement process to clean up the final mesh.
   *
   * @param lcc The Linear Cell Complex containing the mesh
   * @param func Trimming function that determines whether a volume should be kept
   *             (returns true to keep, false to remove)
   */
  void trim_excedent_volumes(LCC& lcc, TrimmingFunction func){
    auto volumes = lcc.one_dart_per_cell<3>();
    for (auto it = volumes.begin(); it != volumes.end(); it++){
      if (func(lcc, it)) continue;
      lcc.remove_cell<3>(it);
    }
  }

  /**
   * @brief Returns the 8 hexahedral volumes (3-cells) surrounding a given node (vertex).
   *
   * This function returns the 8 volumes (hexahedra) that share the specified node in a regular grid.
   * The order of the returned volumes is as follows:
   *   0: left-top-back
   *   1: right-top-back
   *   2: right-top-front
   *   3: left-top-front
   *   4: left-bottom-back
   *   5: right-bottom-back
   *   6: right-bottom-front
   *   7: left-bottom-front
   *
   * The arrangement corresponds to the following dual cube:
   *
   *        z
   *        ^　　　　　y
   *        |  ┐
   *        | /
   *        |/
   *        +------> x
   *
   *      0-------1
   *     /|      /|
   *    3-------2 |
   *    | 4-----|-5
   *    |/      |/
   *    7-------6
   *
   * Each number represents the index in the returned array.
   *
   * @param lcc The Linear Cell Complex
   * @param dart A dart handle representing the node (vertex)
   * @return std::array<Dart_handle, 8> Array of dart handles to the 8 surrounding volumes
   */
  std::array<Dart_handle, 8> volumes_around_node(LCC& lcc, Dart_handle dart) {
    std::array<Dart_handle, 8> volumes = {dart, lcc.beta(dart, 3), lcc.beta(dart, 3, 2, 3), lcc.beta(dart, 2, 3),
                                    lcc.beta(dart, 0, 2, 3), lcc.beta(dart, 3, 1, 2, 3), lcc.beta(dart, 3, 2, 3, 1, 2, 3), lcc.beta(dart, 2, 3, 0, 2, 3)};
    return volumes;
  }

  std::pair<Vector, Vector> get_orthogonal_vectors(const Vector& v) {
    Vector base = (std::abs(v.x()) < std::abs(v.y())) ? Vector(1, 0, 0) : Vector(0, 1, 0);

    Vector v1 = CGAL::cross_product(v, base);
    Vector v2 = CGAL::cross_product(v, v1);

    return {v1/std::sqrt(v1*v1), v2/std::sqrt(v2*v2)};
  }
}




#endif