// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)    Andreas Fabri, Laurent Rineau

#ifndef CGAL_TRIANGULATION_DATA_STRUCTURE_3_FWD_H
#define CGAL_TRIANGULATION_DATA_STRUCTURE_3_FWD_H

#include <CGAL/license/TDS_3.h>

#include <CGAL/tags.h>
#include <CGAL/Triangulation_ds_cell_base_3.h>
#include <CGAL/Triangulation_ds_vertex_base_3.h>

namespace CGAL {

template < class Vb = Triangulation_ds_vertex_base_3<>,
           class Cb = Triangulation_ds_cell_base_3<>,
           class Concurrency_tag_ = Sequential_tag,
           class Storage_policy = Handle_tag
>
class Triangulation_data_structure_3;

} // namespace CGAL

#endif // CGAL_TRIANGULATION_DATA_STRUCTURE_3_FWD_H
