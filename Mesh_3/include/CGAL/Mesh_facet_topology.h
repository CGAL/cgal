// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_MESH_FACET_TOPOLOGY_H
#define CGAL_MESH_FACET_TOPOLOGY_H

#include <CGAL/license/Mesh_3.h>


namespace CGAL {

enum Mesh_facet_topology {
  FACET_VERTICES_ON_SURFACE = 1,
  FACET_VERTICES_ON_SAME_SURFACE_PATCH = 2,
  FACET_VERTICES_ON_SAME_SURFACE_PATCH_WITH_ADJACENCY_CHECK = 3,
  MANIFOLD_WITH_BOUNDARY = 8,
  NO_BOUNDARY = 16,
  MANIFOLD = 24
};

} // end namespace CGAL

#endif // CGAL_MESH_FACET_TOPOLOGY_H
