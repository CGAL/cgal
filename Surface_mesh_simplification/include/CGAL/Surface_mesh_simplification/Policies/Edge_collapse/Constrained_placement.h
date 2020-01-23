// Copyright (c) 2014  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot <sebastien.loriot@cgal.org>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_CONSTRAINED_PLACEMENT_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_CONSTRAINED_PLACEMENT_H

#include <CGAL/license/Surface_mesh_simplification.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>"
#define CGAL_DEPRECATED_MESSAGE_DETAILS "See the Named Parameter `constrain_geometry` for more information."
#include <CGAL/internal/deprecation_warning.h>

#include <CGAL/Surface_mesh_simplification/internal/Placement_wrappers.h>

namespace CGAL {
namespace Surface_mesh_simplification {

template<class BasePlacement, class EdgeIsConstrainedMap>
class Constrained_placement
  : public internal::Constrained_placement<BasePlacement, EdgeIsConstrainedMap>
{
public:
  Constrained_placement(const EdgeIsConstrainedMap map = EdgeIsConstrainedMap(),
                        const BasePlacement& base = BasePlacement())
    : internal::Constrained_placement<BasePlacement, EdgeIsConstrainedMap>(map, base)
  {}
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_CONSTRAINED_PLACEMENT_H
