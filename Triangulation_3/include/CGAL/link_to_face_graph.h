// Copyright (c) 2017  GeometryFactory (France).
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
// Author(s)     : Andreas Fabri
//


#ifndef CGAL_LINK_TO_FACE_GRAPH_H
#define CGAL_LINK_TO_FACE_GRAPH_H

#include <CGAL/license/Triangulation_3.h>


#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
#include <CGAL/array.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/helpers.h>

namespace CGAL {


template<class Triangulation_3,class FG>
typename boost::graph_traits<FG>::vertex_descriptor
link_to_face_graph(const Triangulation_3& t,
                   typename Triangulation_3::Vertex_handle vh,
                   FG& fg, 
                   bool no_infinite_faces = true)
{
  typedef typename Triangulation_3::Cell_handle Cell_handle;
  typedef typename Triangulation_3::Vertex_handle Vertex_handle;
  typedef typename boost::graph_traits<FG>::vertex_descriptor vertex_descriptor;

  clear(fg);
  vertex_descriptor inf;
  vertex_descriptor nullvertex = boost::graph_traits<FG>::null_vertex();
  fg.clear();
  typedef boost::unordered_map<Vertex_handle, vertex_descriptor> Vertex_map;
  Vertex_map vertex_map;
  std::vector<Cell_handle>  cells;
  t.incident_cells(t.infinite_vertex(),std::back_inserter(cells));
  CGAL::cpp11::array<vertex_descriptor,3> face;

  typename boost::property_map<FG, CGAL::vertex_point_t>::type vpm
    = get(CGAL::vertex_point, fg);

  BOOST_FOREACH(Cell_handle ch, cells){
    bool infinite_face = false;
    int vhi = ch->index(vh);
    for(int i=0; i<3; i++){
      int j = Triangulation_3::vertex_triple_index(vhi,i);
      Vertex_handle vhj = ch->vertex(j);
      if(no_infinite_faces && t.is_infinite(vhj)){
        infinite_face = true;
      } else {
        std::pair<typename Vertex_map::iterator,bool> res 
          = vertex_map.insert(std::make_pair(vhj,nullvertex));
        if(res.second){
          res.first->second = add_vertex(fg);
          put(vpm, res.first->second, vhj->point());
          if(t.is_infinite(vhj)){
            inf = res.first->second;
          }
        }
        face[i] = res.first->second;
      }
    }
    if(! infinite_face){
      Euler::add_face(face,fg);
    }
  }
  return inf;
}

} //namespace CGAL

#endif //CGAL_LINK_TO_FACE_GRAPH_H
