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
// Author(s)     : Antonio Gomes, Dmitry Anisimov
//

#ifndef CGAL_BARYCENTRIC_ENUM_3_H
#define CGAL_BARYCENTRIC_ENUM_3_H

// #include <CGAL/license/Barycentric_coordinates_3.h>

namespace CGAL {

/*!
  \ingroup PkgBarycentricCoordinates3Ref

  The namespace `Barycentric_coordinates` contains implementations of all
  generalized barycentric coordinates: 2D, 3D, related enumerations, etc.
*/
namespace Barycentric_coordinates {

enum class Computation_policy_3 {

  FAST = 0,
  FAST_WITH_EDGE_CASES = 1

};


} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_ENUM_3_H
