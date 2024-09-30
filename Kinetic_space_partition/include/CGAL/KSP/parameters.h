// Copyright (c) 2021 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot, Dmitry Anisimov

#ifndef CGAL_KSP_PARAMETERS_H
#define CGAL_KSP_PARAMETERS_H

#include <CGAL/license/Kinetic_space_partition.h>

namespace CGAL {
namespace KSP {
namespace internal {

template<typename FT>
struct Parameters_3 {

  unsigned int k = 1; // k intersections

  // Octree parameters for subdivison of input data into kinetic subpartitions
  unsigned int max_octree_depth = 3;
  unsigned int max_octree_node_size = 40;

  FT bbox_dilation_ratio = FT(11) / FT(10); // ratio to enlarge bbox

  bool reorient_bbox = false; // true - optimal bounding box, false - axis aligned

  // All files are saved in the current build directory.
  bool verbose = false; // print basic verbose information
  bool debug = false; // print all steps and substeps + export initial and final configurations

  // See also global tolerance inside utils.h! (set to 0)
  Parameters_3(const bool v = true, const bool d = false) :
    verbose(v), debug(d) { }
};

} // namespace internal
} // namespace KSP
} // namespace CGAL

#endif // CGAL_KSP_PARAMETERS_H
