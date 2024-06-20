// Copyright (c) 2019-2024  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_CELL_BASE_3_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_CELL_BASE_3_H

#include <CGAL/Triangulation_cell_base_3.h>

namespace CGAL {

/**
 * @ingroup PkgCT_3Classes
 * @brief Cell base class for the 3D Constrained Delaunay Triangulation.
 *
 * This class is derived from the `Triangulation_cell_base_3` class and provides additional functionality
 * required by `make_constrained_Delaunay_triangulation_3()`.
 *
 * @tparam Gt The geometric traits class, which must be a model of `ConstrainedDelaunayTriangulationTraits_3`.
 *         It should be the same as the geometric traits class of the triangulation.
 * @tparam Cb The base class for the cell, which must be a model of `TriangulationCellBase_3`.
 *
 * @cgalModels{ConstrainedDelaunayTriangulationCellBase_3, SimplicialMeshCellBase_3, RemeshingCellBase_3}
 *
 * \note This cell base class also models the `SimplicialMeshCellBase_3` and `RemeshingCellBase_3` concepts, allowing the use of functionality from \ref Chapter_Tetrahedral_Remeshing "Tetrahedral Remeshing" and \ref Chapter_3D_Simplicial_Mesh_Data_Structure "3D Simplicial Mesh Data Structures", if the corresponding vertex base also models the right concepts.
 * \todo After discussion with Jane. The note above is wrong. There should be a second pair of Vb/Cb, designed to model the concepts of simplicial mesh and remeshing.
 *
 * \sa `CGAL::Constrained_Delaunay_triangulation_vertex_base_3`
 */
template < typename Gt, typename Cb = Triangulation_cell_base_3<Gt> >
class Constrained_Delaunay_triangulation_cell_base_3 : public Cb {
public:
};

} // namespace CGAL

#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_CELL_BASE_3_H
