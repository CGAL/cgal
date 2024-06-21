// Copyright (c) 2024  GeometryFactory Sarl (France).
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

#ifndef CGAL_CT_3_TYPES_H
#define CGAL_CT_3_TYPES_H

#include <CGAL/license/Triangulation_3.h>

#include <CGAL/config.h>

namespace CGAL {

/**
 * @addtogroup PkgCT_3Classes
 * @typedef CDT_3_face_index
 * Integral type to store the index of constraints.
 * @see `Constrained_Delaunay_triangulation_cell_data_3`
 * @see `Constrained_Delaunay_triangulation_vertex_base_3`
 *
 */
using CDT_3_face_index = int; // must be signed

} // namespace CGAL

#endif // CGAL_CT_3_TYPES_H
