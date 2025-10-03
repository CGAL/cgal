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
#ifndef CGAL_HEXMESHING_ADJACENT_FACE_ON_PLANE_H
#define CGAL_HEXMESHING_ADJACENT_FACE_ON_PLANE_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_grid.h>
#include <boost/container/static_vector.hpp>

namespace CGAL::internal::Hexmeshing {
  /**
   * @brief Internal function to find the adjacent face on a given plane
   *
   * This function searches for a face adjacent to the given edge that lies on the specified plane.
   * It traverses through the mesh using beta operations to find faces that share the edge
   * and checks if they belong to the target plane. The function may fail if the adjacent face
   * is not directly accessible from beta3 to the current face. To complete the search, you should
   * call this function with both edge and lcc.beta<3>(edge)
   *
   * @param lcc The Linear Cell Complex
   * @param plane The plane normal direction to search for
   * @param edge A dart handle representing the edge to find adjacent faces for
   * @return A dart handle to the adjacent face on the plane, or null if not found
   */
  Dart_handle __adjacent_face_on_plane(LCC& lcc, PlaneNormal plane, Dart_handle edge){
    Dart_handle other_face = lcc.beta(edge, 2);
    auto this_face_handle = lcc.attribute<2>(edge);
    auto other_face_handle = lcc.attribute<2>(other_face);

    assert(this_face_handle != nullptr);

    // Early exit if the beta(2) face is already on the plane
    if (other_face_handle != nullptr
      && other_face_handle->info().plane[plane])
      return other_face;

    bool found = false;
    for (int i = 0; i < 5 && !found; i++){

      if (lcc.is_free<3>(other_face)) {
        // We hit the grid border, there are no adjacent face possible to that edge
        break;
      }

      other_face = lcc.beta(other_face, 3, 2);
      other_face_handle = lcc.attribute<2>(other_face);
      auto other_vol_handle = lcc.attribute<3>(other_face);

      // Exit if we fall back on the same face or outside of the domain
      if (this_face_handle == other_face_handle)
        return lcc.null_dart_descriptor;

      found = other_face_handle != nullptr
        && other_face_handle->info().plane[plane];
    }

    return found ? other_face : lcc.null_dart_descriptor;
  }

  /**
   * @brief Internal function to find the adjacent face on a given plane with additional volume collection
   *
   * This overload of __adjacent_face_on_plane additionally collects volumes that couldn't be reached
   * directly from the plane during the search. It searches for a face adjacent to the given edge
   * that lies on the specified plane, and stores any encountered volumes that are inside the
   * refinement domain but not directly accessible from the plane.
   *
   * @param lcc The Linear Cell Complex
   * @param plane The plane normal direction to search for
   * @param edge A dart handle representing the edge to find adjacent faces for
   * @param additional_volumes Vector to store volumes encountered during the search
   * @return A dart handle to the adjacent face on the plane, or null if not found
   */
  Dart_handle __adjacent_face_on_plane(LCC& lcc, PlaneNormal plane, Dart_handle edge, std::vector<Dart_handle>& additional_volumes){
    Dart_handle other_face = lcc.beta(edge, 2);
    auto this_face_handle = lcc.attribute<2>(edge);
    auto other_face_handle = lcc.attribute<2>(other_face);
    boost::container::static_vector<Dart_handle, 5> __additional_volumes;

    assert(this_face_handle != nullptr);

    // Early exit if the beta(2) face is already on the plane
    if (other_face_handle != nullptr && other_face_handle->info().plane[plane])
      return other_face;

    bool found = false;
    for (int i = 0; i < 5 && !found; i++){

      if (lcc.is_free<3>(other_face)) {
        // We hit the grid border, there are no adjacent face possible to that edge
        break;
      }

      other_face = lcc.beta(other_face, 3, 2);
      other_face_handle = lcc.attribute<2>(other_face);
      auto other_vol_handle = lcc.attribute<3>(other_face);

      // Exit if we fall back on the same face or outside of the domain
      if (this_face_handle == other_face_handle)
        break;

      found = other_face_handle != nullptr && other_face_handle->info().plane[plane];

      // Add the incident volume only if it is inside the refinement domain
      if (!found && other_vol_handle != nullptr)
        __additional_volumes.push_back(other_face);
    }

    additional_volumes.insert(additional_volumes.end(), __additional_volumes.begin(), __additional_volumes.end());

    return found ? other_face : lcc.null_dart_descriptor;
  }

  /**
   * @brief Finds the adjacent face on a given plane
   *
   * This function searches for a face adjacent to the given edge that lies on the specified plane.
   * It uses the internal __adjacent_face_on_plane function and handles cases where the adjacent face
   * might not be directly accessible from the current orientation. If the first search fails,
   * it tries again with the opposite orientation (beta<3> of the edge).
   *
   * @param lcc The Linear Cell Complex
   * @param plane The plane normal direction to search for
   * @param edge A dart handle representing the edge to find adjacent faces for
   * @return A dart handle to the adjacent face on the plane, or null if not found
   */
  Dart_handle adjacent_face_on_plane(LCC& lcc, PlaneNormal plane, Dart_handle edge){
    Dart_handle other_face = __adjacent_face_on_plane(lcc, plane, edge);

    // This is needed to iterate on 3 templates that are on a grid border, missing a connected volume
    // And preventing iteration the first way
    if (other_face == lcc.null_dart_descriptor && !lcc.is_free<3>(edge))
      other_face = __adjacent_face_on_plane(lcc, plane, lcc.beta<3>(edge));

    return other_face;
  }

  /**
   * @brief Finds the adjacent face on a given plane with additional volume collection
   *
   * This overload of adjacent_face_on_plane additionally collects volumes that couldn't be reached
   * directly from the plane during the search. It uses the internal __adjacent_face_on_plane function
   * with volume collection and handles cases where the adjacent face might not be directly accessible
   * from the current orientation. If the first search fails, it tries again with the opposite orientation.
   *
   * @param lcc The Linear Cell Complex
   * @param plane The plane normal direction to search for
   * @param edge A dart handle representing the edge to find adjacent faces for
   * @param additionnal_volumes Pointer to vector to store volumes encountered during the search
   * @return A dart handle to the adjacent face on the plane, or null if not found
   */
  Dart_handle adjacent_face_on_plane(LCC& lcc, PlaneNormal plane, Dart_handle edge, std::vector<Dart_handle>* additionnal_volumes){
    Dart_handle other_face = __adjacent_face_on_plane(lcc, plane, edge, *additionnal_volumes);

    if (other_face == lcc.null_dart_descriptor && !lcc.is_free<3>(edge))
      other_face = __adjacent_face_on_plane(lcc, plane, lcc.beta<3>(edge), *additionnal_volumes);

    return other_face;
  }
}

#endif // CGAL_HEXMESHING_ADJACENT_FACE_ON_PLANE_H