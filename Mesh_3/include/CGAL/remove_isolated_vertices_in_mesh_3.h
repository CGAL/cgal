// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois
//
//******************************************************************************
// File Description : remove_isolated_vertices_in_mesh_3 function definition.
//******************************************************************************

#ifndef CGAL_REMOVE_ISOLATED_VERTICES_IN_MESH_3_H
#define CGAL_REMOVE_ISOLATED_VERTICES_IN_MESH_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

namespace CGAL {

/*!
  \ingroup PkgMesh3Functions
  The tetrahedral mesh generation algorithm implemented in `CGAL::make_mesh_3()`
  and `CGAL::refine_mesh_3()` does not guarantee that all the points inserted
  by the algorithm are actually present in the final mesh.
  In most cases, all points are used, but if the geometry of the object
  has small features, compared to the size of the simplices (triangles and tetrahedra),
  it might be that the Delaunay facets that are selected in the restricted Delaunay
  triangulation miss some vertices of the triangulation.
  This function removes these so-called "isolated" vertices, that belong to the
  triangulation but not to any cell of the `C3T3`, from the triangulation.

  \tparam C3T3 is required to be a model of the concept `MeshComplex_3InTriangulation_3`.
  The argument `c3t3` is passed by
  reference as its underlying triangulation is modified by the function.
*/
template <typename C3T3>
void
remove_isolated_vertices_in_mesh_3(C3T3& c3t3)
{
  c3t3.remove_isolated_vertices();
}


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_REMOVE_ISOLATED_VERTICES_IN_MESH_3_H
