// Copyright (c) 2016 GeometryFactory (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_PMP_INTERNAL_COREFINEMENT_INTERSECTION_TYPE_H
#define CGAL_PMP_INTERNAL_COREFINEMENT_INTERSECTION_TYPE_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <boost/graph/graph_traits.hpp>
#include <CGAL/property_map.h>

namespace CGAL{
namespace Corefinement{

enum Intersection_type {ON_FACE,ON_EDGE,ON_VERTEX,EMPTY,COPLANAR_TRIANGLES};

template <class TriangleMesh, class Exact_kernel>
struct Coplanar_intersection{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  Intersection_type type_1,type_2;   //intersection type for 1st and 2nd faces
  halfedge_descriptor info_1,info_2; //halfedge providing primitive indicated by type_1 and type_2
  typename Exact_kernel::Point_3 point;


  Coplanar_intersection()
  : type_1(EMPTY), type_2(EMPTY)
  , info_1(GT::null_halfedge()), info_2(info_1)
  {}

};

} } // CGAL::Corefinement

#endif // CGAL_PMP_INTERNAL_COREFINEMENT_INTERSECTION_TYPE_H
