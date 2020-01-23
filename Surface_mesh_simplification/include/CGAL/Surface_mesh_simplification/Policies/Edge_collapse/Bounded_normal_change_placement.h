// Copyright (c) 2017  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_NORMAL_CHANGE_PLACEMENT_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_NORMAL_CHANGE_PLACEMENT_H

#include <CGAL/license/Surface_mesh_simplification.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>"
#define CGAL_DEPRECATED_MESSAGE_DETAILS "See the Named Parameter `max_normal_angle_change` for more information."

#include <CGAL/internal/deprecation_warning.h>
#include <CGAL/Surface_mesh_simplification/internal/Placement_wrappers.h>

namespace CGAL {
namespace Surface_mesh_simplification {

template<class GetPlacement>
class Bounded_normal_change_placement
    : public internal::Bounded_normal_change_placement<GetPlacement>
{
public:
   Bounded_normal_change_placement(const GetPlacement& get_placement = GetPlacement())
     :internal::Bounded_normal_change_placement<GetPlacement>(CGAL_PI/4.0, get_placement)
   {}
};
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_NORMAL_CHANGE_PLACEMENT_H
