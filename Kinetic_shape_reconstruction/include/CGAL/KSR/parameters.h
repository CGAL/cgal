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

#ifndef CGAL_KSR_PARAMETERS_H
#define CGAL_KSR_PARAMETERS_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

namespace CGAL {
namespace KSR {

template<typename FT>
struct Parameters_3 {

  unsigned int k = 1; // k intersections
  unsigned int n = 0; // n subdivisions, not implemented yet

  FT bbox_dilation_ratio = FT(11) / FT(10); // ratio to enlarge bbox
  FT angle_tolerance = FT(5);
  FT distance_tolerance =  FT(5) / FT(10); // distance tolerance between planes

  bool reorient_bbox = false; // true - optimal bounding box, false - axis aligned

  // All files are saved in the current build directory.
  bool verbose    =  true; // print basic verbose information
  bool debug      = false; // print all steps and substeps + export initial and final configurations

  // See also global tolerance inside utils.h!
  Parameters_3(const bool v = true, const bool d = false) :
  verbose(v), debug(d) { }
};

} // namespace KSR
} // namespace CGAL

#endif // CGAL_KSR_PARAMETERS_H
