// Copyright (c) 2011  GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_CONVEX_HULL_3_TO_FACE_GRAPH_3_H
#define CGAL_CONVEX_HULL_3_TO_FACE_GRAPH_3_H

#include <CGAL/license/Convex_hull_3.h>


#include <CGAL/link_to_face_graph.h>

namespace CGAL {


template<class Triangulation_3,class PolygonMesh>
void convex_hull_3_to_face_graph(const Triangulation_3& T,PolygonMesh& P){
  link_to_face_graph(T,T.infinite_vertex(), P);
}

} //namespace CGAL

#endif //CGAL_CONVEX_HULL_3_TO_FACE_GRAPH_3_H
