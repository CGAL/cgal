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

#ifndef CGAL_TET_ADAPTIVE_REMESHING_CELL_BASE_3_H
#define CGAL_TET_ADAPTIVE_REMESHING_CELL_BASE_3_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Simplicial_mesh_cell_base_3.h>

namespace CGAL
{
namespace Tetrahedral_remeshing
{

/*!
\ingroup PkgTetrahedralRemeshingClasses

The class `Remeshing_cell_base_3` is a model of the concept `MeshCellBase_3`.
It is designed to serve as cell base class for the 3D triangulation
used in the tetrahedral remeshing process.

\todo indices

\cgalModels `SimplicialMeshCellBase_3`

*/
template<typename Subdomain_index = int,
         typename Surface_patch_index = int,
         typename TDS = void>
using Remeshing_cell_base_3
  = CGAL::Simplicial_mesh_cell_base_3<Subdomain_index,
                                      Surface_patch_index,
                                      TDS>;

}//end namespace Tetrahedral_remeshing
}//end namespace CGAL

#endif //CGAL_TET_ADAPTIVE_REMESHING_CELL_BASE_3_H
