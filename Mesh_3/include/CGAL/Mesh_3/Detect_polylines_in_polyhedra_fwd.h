// Copyright (c) 2010 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau
//

#ifndef CGAL_DETECT_POLYLINES_IN_POLYHEDRA_FWD_H
#define CGAL_DETECT_POLYLINES_IN_POLYHEDRA_FWD_H

#include <CGAL/license/Mesh_3.h>


namespace CGAL { namespace Mesh_3 {

template <typename Polyhedron>
struct Detect_polylines;

template <typename Polyhedron,
          typename Polyline_and_context,
          typename Polylines_output_iterator>
Polylines_output_iterator
detect_polylines(Polyhedron* pMesh,
                 Polylines_output_iterator out_it);

} // end namespace CGAL::Mesh_3
} // end namespace CGAL


#endif // CGAL_DETECT_POLYLINES_IN_POLYHEDRA_FWD_H
