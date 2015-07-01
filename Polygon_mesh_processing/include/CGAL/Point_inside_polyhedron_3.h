// Copyright (c) 2013,2014,2015 GeometryFactory (France).
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
// Author(s)     : Jane Tournois

#ifndef CGAL_POINT_INSIDE_POLYHEDRON_H
#define CGAL_POINT_INSIDE_POLYHEDRON_H

#include <CGAL/Side_of_triangle_mesh.h>

namespace CGAL {

// this file is there to make Mean_curvature_flow_skeltonization work 
// with and without PMP
template <class Polyhedron, 
          class Kernel> 
class Point_inside_polyhedron_3 : public  Side_of_triangle_mesh<Polyhedron,Kernel>
{
  typedef Side_of_triangle_mesh<Polyhedron,Kernel> Base;
public:
  Point_inside_polyhedron_3(const Polyhedron& polyhedron)
    : Base(polyhedron)
  {}
};

}

#endif //CGAL_POINT_INSIDE_POLYHEDRON_H 
