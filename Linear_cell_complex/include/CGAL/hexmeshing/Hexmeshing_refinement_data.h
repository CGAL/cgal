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
#ifndef HEXMESHING_REFINEMENT_DATA_H
#define HEXMESHING_REFINEMENT_DATA_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_plane_normal.h>
#include <vector>

namespace CGAL::internal::Hexmeshing {
  /**
   * @brief Data structure for managing mesh refinement operations
   * 
   * This structure contains the necessary data for performing mesh refinement,
   * including iteration direction, marked elements, and collections of elements
   * that need to be refined.
   */
  struct RefinementData {
    PlaneNormal iteration;  ///< Current iteration direction (X, Y, Z) for refinement

    std::vector<Dart_handle> marked_nodes;  ///< Collection of nodes that have been marked during refinement
    std::vector<Dart_handle> faces_of_plane;  ///< Faces lying on the current refinement plane
    std::vector<Dart_handle> additionnal_volumes_found;  ///< Additional volumes discovered during refinement

    std::vector<Dart_handle> volumes_to_refine;  ///< Collection of volumes that need to be refined
    std::vector<Dart_handle> faces_to_refine;  ///< Collection of faces that need to be refined
    std::vector<Dart_handle> partial_templates_to_refine;  ///< Collection of partial templates that need refinement
  };
}



#endif