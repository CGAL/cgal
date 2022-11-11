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

/*!
\ingroup PkgMesh3Enum

The enum `Mesh_facet_topology` is designed to tell which constraints have to
be checked on each surface facet during the mesh refinement process.

\sa `CGAL::Mesh_criteria_3<Tr>`,
\sa `CGAL::Mesh_facet_criteria_3<Tr>`.
*/
enum Mesh_facet_topology
{
  FACET_VERTICES_ON_SURFACE = 1,//!< Each vertex of the facet has
                                //!< to be on the surface, on a curve, or on a corner.
  FACET_VERTICES_ON_SAME_SURFACE_PATCH = 2, //!< The three vertices of a facet belonging
                                            //!< to a surface patch `s` have to be on
                                            //!< the same surface patch `s`, on a curve or on a corner.
  /*!
    The three vertices of a facet belonging to a surface patch `s`
    have to be on the same surface patch `s`, or on a curve
    incident to the surface patch `s` or on a corner incident to the
    surface patch `s`.
  */
  FACET_VERTICES_ON_SAME_SURFACE_PATCH_WITH_ADJACENCY_CHECK = 3
#ifndef DOXYGEN_RUNNING
  ,
  MANIFOLD_WITH_BOUNDARY = 8,
  NO_BOUNDARY = 16,
  MANIFOLD = 24
#endif
};

} // end namespace CGAL

#endif // CGAL_MESH_FACET_TOPOLOGY_H
