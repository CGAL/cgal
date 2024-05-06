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

#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_CELL_DATA_3_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_CELL_DATA_3_H

#include <CGAL/license/Constrained_triangulation_3.h>

#include <CGAL/Constrained_triangulation_3/internal/config.h>

namespace CGAL {
#ifdef DOXYGEN_RUNNING
/*!
 * @brief Internal per-cell data for \cgal 3D constrained Delaunay triangulations
 *
 * This class is an internal detail of the implementation of \cgal 3D constrained Delaunay triangulations.
 *
 * Any model of the `ConstrainedDelaunayTriangulationCellBase_3` concept must include one object of this type
 * as a non-static data member.
 */
struct Constrained_Delaunay_triangulation_cell_data_3 {};
#else // DOXYGEN_RUNNING
#endif // DOXYGEN_RUNNING
} // namespace CGAL

#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_CELL_DATA_3_H
