// Copyright (c) 2011  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// 
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_CONVEX_HULL_3_TO_POLYHEDRON_3_H
#define CGAL_CONVEX_HULL_3_TO_POLYHEDRON_3_H

#include <CGAL/license/Convex_hull_3.h>


#define CGAL_DEPRECATED_HEADER "<CGAL/convex_hull_3_to_polyhedron_3.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/convex_hull_3_to_face_graph.h>"
#include <CGAL/internal/deprecation_warning.h>


#include <CGAL/Polyhedron_3_fwd.h>

namespace CGAL {
  
template<class Triangulation_3,class Polyhedron_3>
CGAL_DEPRECATED void convex_hull_3_to_polyhedron_3(const Triangulation_3& T,Polyhedron_3& P){
  clear(P);
  link_to_face_graph(T,T.infinite_vertex(), P);
}

} //namespace CGAL

#endif //CGAL_CONVEX_HULL_3_TO_POLYHEDRON_3_H
