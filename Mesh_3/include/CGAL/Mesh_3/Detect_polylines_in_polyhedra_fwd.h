// Copyright (c) 2010 GeometryFactory Sarl (France).
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
//
//
// Author(s)     : Laurent Rineau
//

#ifndef CGAL_DETECT_POLYLINES_IN_POLYHEDRA_FWD_H
#define CGAL_DETECT_POLYLINES_IN_POLYHEDRA_FWD_H

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
