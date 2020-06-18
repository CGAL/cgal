// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_INTERNAL_LINDSTROMTURK_PARAMS_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_INTERNAL_LINDSTROMTURK_PARAMS_H

#include <CGAL/license/Surface_mesh_simplification.h>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

struct LindstromTurk_params
{
  LindstromTurk_params()
    :
      m_volume_weight(0.5),
      m_boundary_weight(0.5),
      m_shape_weight(0)
  {}

  LindstromTurk_params(double volume_weight, double boundary_weight, double shape_weight)
    :
      m_volume_weight(volume_weight),
      m_boundary_weight(boundary_weight),
      m_shape_weight(shape_weight)
  {}

  double m_volume_weight;
  double m_boundary_weight;
  double m_shape_weight;
};

} // namespace internal
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_INTERNAL_LINDSTROMTURK_PARAMS_H
