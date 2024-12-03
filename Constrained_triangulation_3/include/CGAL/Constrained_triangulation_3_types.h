// Copyright (c) 2024  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
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
 * @ingroup PkgCT_3Classes
 * \brief Integral type to store the index of constraints.
 * @see `Conforming_constrained_Delaunay_triangulation_cell_data_3`
 * @see `Conforming_constrained_Delaunay_triangulation_vertex_base_3`
 *
 */
using CDT_3_face_index = int; // must be signed

} // namespace CGAL

#endif // CGAL_CT_3_TYPES_H
