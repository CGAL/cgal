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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_COUNT_STOP_PREDICATE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_COUNT_STOP_PREDICATE_H

#include <CGAL/license/Surface_mesh_simplification.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>"
#include <CGAL/Installation/internal/deprecation_warning.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>

#ifndef CGAL_NO_DEPRECATED_CODE

namespace CGAL {
namespace Surface_mesh_simplification {

// Stops when the number of edges left falls below a given number.
template<class TM_>
using Count_stop_predicate CGAL_DEPRECATED = Edge_count_stop_predicate<TM_>;

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_NO_DEPRECATED_CODE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_COUNT_STOP_PREDICATE_H
