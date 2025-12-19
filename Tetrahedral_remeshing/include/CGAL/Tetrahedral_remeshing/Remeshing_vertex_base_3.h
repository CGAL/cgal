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

namespace CGAL {
namespace Tetrahedral_remeshing {

/*!
\ingroup PkgTetrahedralRemeshingClasses

The class `Remeshing_vertex_base_3` is a model of the concept `RemeshingVertexBase_3`.
It is designed to serve as vertex base class for the 3D triangulation
used in the tetrahedral remeshing process.

\tparam Gt is the geometric traits class.
It must be a model of the concept `RemeshingTriangulationTraits_3`.

\tparam Vb is a vertex base class from which `Remeshing_vertex_base_3` derives.
It must be a model of the concept `SimplicialMeshVertexBase_3`.

\cgalModels{RemeshingVertexBase_3,SimplicialMeshVertexBase_3}
*/
template<typename Gt,
         typename Vb = CGAL::Simplicial_mesh_vertex_base_3<Gt,
                         int /*Subdomain_index*/,
                         int /*Surface_patch_index*/,
                         int /*Curve_index*/,
                         int /*Corner_index*/> >
class Remeshing_vertex_base_3
  : public Vb
{
public:
  template <typename TDS2>
  struct Rebind_TDS
  {
    using Vb2 = typename Vb::template Rebind_TDS<TDS2>::Other;
    using Other = Remeshing_vertex_base_3<Gt, Vb2>;
  };

public:
  using Vb::Vb; // constructors
};

} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif //CGAL_TET_ADAPTIVE_REMESHING_VERTEX_BASE_3_H
