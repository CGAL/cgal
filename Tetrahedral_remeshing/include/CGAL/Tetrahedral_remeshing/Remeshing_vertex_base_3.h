// Copyright (c) 2020 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Noura Faraj, Jean-Marc Thiery, Tamy Boubekeur

#ifndef CGAL_TET_ADAPTIVE_REMESHING_VERTEX_BASE_3_H
#define CGAL_TET_ADAPTIVE_REMESHING_VERTEX_BASE_3_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Simplicial_mesh_vertex_base_3.h>

namespace CGAL
{
namespace Tetrahedral_remeshing
{

/*!
\ingroup PkgTetrahedralRemeshingClasses

The class `Remeshing_vertex_base_3` is a model of the concept `MeshVertexBase_3`.
It is designed to serve as vertex base class for the 3D triangulation
used in the tetrahedral remeshing process.

\tparam Gt is the geometric traits class.
It has to be a model of the concept `RemeshingTriangulationTraits_3`.

\tparam Vb is a vertex base class from which `Remeshing_vertex_base_3` derives.
It must be a model of the `TriangulationVertexBase_3` concept.
It has the default value `Triangulation_vertex_base_3<Gt>`.

\cgalModels{RemeshingVertexBase_3,SimplicialMeshVertexBase_3}
*/
template<typename GT,
  typename Subdomain_index = int,
  typename Surface_patch_index = int,
  typename Curve_index = int,
  typename Corner_index = int,
  typename Vb = CGAL::Triangulation_vertex_base_3<GT> >
using Remeshing_vertex_base_3
  = CGAL::Simplicial_mesh_vertex_base_3<GT,
                                        Subdomain_index,
                                        Surface_patch_index,
                                        Curve_index,
                                        Corner_index,
                                        Vb>;

}//end namespace Tetrahedral_remeshing

}//end namespace CGAL

#endif //CGAL_TET_ADAPTIVE_REMESHING_VERTEX_BASE_3_H
