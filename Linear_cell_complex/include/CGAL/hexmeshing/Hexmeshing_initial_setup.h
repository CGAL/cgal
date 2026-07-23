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
#ifndef HEXMESHING_INITIAL_SETUP_H
#define HEXMESHING_INITIAL_SETUP_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>

namespace CGAL::internal::Hexmeshing {
  /**
   * @brief Sets up the initial plane structure for the hexahedral mesh refinement algorithm
   *
   * This function initializes the plane structure that will be used throughout the refinement
   * process. It creates and organizes faces into planes along the X, Y, and Z axes, establishing
   * the foundation for the refinement algorithm's plane-based approach.
   *
   * The function performs two main operations:
   *
   * 1. **Plane Extraction**: For each axis (X, Y, Z), it extracts the first faces of each plane
   *    by traversing the mesh structure using beta operations. It uses helper functions:
   *    - `__first_plane`: Finds the starting face for each plane direction
   *    - `__next_plane`: Navigates to the next plane in the sequence
   *    The extraction iterates through all planes in each direction based on the grid dimensions.
   *
   * 2. **Attribute Creation**: For all faces in all planes, it creates and initializes face
   *    attributes with the following properties:
   *    - Sets the appropriate plane bit in the plane bitset
   *    - Assigns a plane_id based on the plane's position in the sequence
   *    - Sets cc_id (connected component ID) to 0 for all faces
   *
   * The function uses `plane_for_each_face` to traverse all faces in each plane and ensure
   * proper attribute initialization. This setup is crucial for the refinement algorithm to
   * properly identify and process faces during the refinement stages.
   *
   * Note: The last odd plane in each direction is not considered in the current implementation,
   * as indicated by the commented code at the end of the function.
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and grid information
   */
  template <typename HexData>
  void setup_initial_planes(HexData& hdata){
    LCC& lcc = hdata.lcc;

    auto __first_plane = [](LCC& lcc, PlaneNormal plane){
      switch (plane){
        case X: return lcc.beta(lcc.first_dart(), 0, 2);
        case Y: return lcc.beta(lcc.first_dart(), 2, 1);
        case Z: return lcc.first_dart();
        default:;
      }

      CGAL_assertion_msg(false, "Unexpected value for plane");
      return lcc.null_dart_descriptor;
    };

    auto __next_plane = [](LCC& lcc, Dart_handle start_plane, PlaneNormal plane){
      return lcc.beta(start_plane, 1, 2, 3, 2, 1);
    };

    // Extract the first faces of each plane
    // Note: The last odd plane is not considered
    for (int p = 0; p < 3; p++){
      PlaneNormal plane = (PlaneNormal)p;

      auto& plane_set = hdata.first_face_of_planes[p];

      Dart_handle start_plane = __first_plane(hdata.lcc, plane);

      for (int z = 0; z < hdata.grid.dims[p]; z++){
        std::vector<Dart_handle> starts;
        starts.push_back(lcc.beta(start_plane, 0, 2));

        start_plane = __next_plane(lcc, start_plane, plane);

        plane_set.push_back(std::move(starts));
      }
    }

    // Create attributes of all faces of all planes
    // To be able to iterate properly at later stages
    for (int p = 0; p < 3; p++){
      auto& plane_set = hdata.first_face_of_planes[p];
      for (int i = 0; i < plane_set.size(); i++) {
        plane_for_each_face(lcc, plane_set[i],
          [&](Dart_handle face, auto& edges){
            auto& face_attr = get_or_create_attr<2>(lcc, face)->info();
            face_attr.plane[p] = true;
            face_attr.plane_id = i;
            face_attr.cc_id = 0;
          },
          [&](Dart_handle edge){
            return lcc.beta(edge, 2, 3, 2);
          });
      }

       // Last odd plane
       // int size = hdata.first_face_of_planes[p].size();
       // Dart_handle start = hdata.first_face_of_planes[p][size-1][0];
       // start = lcc.beta(start, 2, 1, 1, 2);

       // plane_for_each_face(lcc, start,
       //  [&](Dart_handle face){
       //    auto& face_attr = get_or_create_attr<2>(lcc, face)->info();
       //    face_attr.plane[p] = true;
       //    face_attr.cc_id = 0;
       //  },
       //  [&](Dart_handle edge){
       //    return lcc.beta(edge, 2, 3, 2);
       //  });
    }
  }

  /**
   * @brief Performs the initial setup for the hexahedral mesh refinement algorithm
   *
   * This function performs the essential initialization steps required before starting
   * the refinement process. It establishes the foundational structure that the refinement
   * algorithm will use throughout its execution.
   *
   * The function performs two main operations:
   *
   * 1. **Plane Setup**: Calls `setup_initial_planes` to initialize the plane structure
   *    along the X, Y, and Z axes. This creates the organizational framework for
   *    the refinement algorithm's plane-based approach.
   *
   * 2. **Cell Identification**: Iterates through all volumes in the Linear Cell Complex
   *    and identifies which cells should be marked for refinement. For each volume:
   *    - Creates a volume attribute if it doesn't exist
   *    - Uses the provided `cellIdentifier` function to determine if the volume
   *      should be marked as identified for refinement
   *    - Sets the volume type to `VolumeType::IDENTIFIED` if the cell should be refined
   *
   * This initialization is crucial because it establishes which cells will be the
   * starting points for the refinement process. The identified cells will be used
   * to determine the refinement patterns and guide the subsequent refinement operations.
   *
   * @tparam HexData Type of the hexahedral meshing data structure
   * @param hdata Hexahedral meshing data containing the Linear Cell Complex and grid information
   * @param cellIdentifier Function that determines whether a cell should be identified
   *                       for refinement based on its position and properties
   */
  template <typename HexData>
  void initial_setup(HexData& hdata, MarkingFunction& cellIdentifier){
    LCC& lcc = hdata.lcc;
    setup_initial_planes(hdata);

    // Mark initial identified cells
    auto volumes = lcc.one_dart_per_cell<3>();
    for (auto dart = volumes.begin(), end = volumes.end(); dart != end; dart++){
      // Create a 3-attr for all 3-cells in the LCC
      auto& vol_attr = get_or_create_attr<3>(lcc, dart)->info();

      // Mark those who are identified
      if (cellIdentifier(lcc, dart))
        vol_attr.type = VolumeType::IDENTIFIED;
    }
  }
}




#endif