// Copyright (c) 2013  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Xiang Gao <gaox@ethz.ch>
//

#ifndef CGAL_MCFSKEL_TRACK_CORRESPONDENCE_H
#define CGAL_MCFSKEL_TRACK_CORRESPONDENCE_H

/// @cond CGAL_DOCUMENT_INTERNAL

/** 
 * @file Track_correspondence_visitor.h
 * @brief This file contains the visitor class to track the correspondent vertices
 * during edge collapse.
 *
 * The visitor class track vertices that are collapsed onto the given vertex.
 * It is needed when using simplification package to do the edge collapse.
 * 
 */

// Visitor base
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
#include <cmath>

namespace SMS = CGAL::Surface_mesh_simplification;

namespace CGAL {
namespace internal {

template<class TriangleMesh, class TriangleMeshPointPMap>
struct Track_correspondence_visitor : SMS::Edge_collapse_visitor_base<TriangleMesh>
{
  // TriangleMesh types
  typedef typename TriangleMesh::Traits         Kernel;
  typedef typename Kernel::Point_3               Point;
  typedef typename TriangleMesh::Vertex_handle  Vertex_handle;

  // Mesh simplification types
  typedef SMS::Edge_profile<TriangleMesh>       Profile;

  Track_correspondence_visitor(){}

  Track_correspondence_visitor(TriangleMeshPointPMap* point_pmap,
                               std::map<int, std::vector<int> >* corr,
                               int max_id) :
    hg_point_pmap(point_pmap),
    corr(corr), 
    max_id(max_id), 
    is_medially_centered(false)
  {}

  Track_correspondence_visitor(TriangleMeshPointPMap* point_pmap,
                       std::map<int, std::vector<int> >* corr,
                       std::map<int, int>* poles,
                       std::vector<Point>* cell_dual,
                       int max_id) :
    hg_point_pmap(point_pmap),
    corr(corr), 
    max_id(max_id), 
    is_medially_centered(true),
    poles(poles), 
    cell_dual(cell_dual)
    {}

  // Called AFTER each edge has been collapsed
  void OnCollapsed(Profile const& edge, Vertex_handle v)
  {
    Vertex_handle v0 = edge.v0();
    Vertex_handle v1 = edge.v1();
    int id0 = v0->id();
    int id1 = v1->id();
    int vid = v->id();
    int from, to;
    if (id0 == vid)
    {
      from = id1;
      to = id0;
    }
    else if (id1 == vid)
    {
      from = id0;
      to = id1;
    }

    if ((*corr).find(to) == (*corr).end())
    {
      (*corr)[to] = std::vector<int>();
    }
    // only track vertex in the original mesh
    if (from < max_id)
    {
      (*corr)[to].push_back(from);
    }
    std::map<int, std::vector<int> >::iterator iter = (*corr).find(from);
    if (iter != (*corr).end())
    {
      for (size_t i = 0; i < (iter->second).size(); ++i)
      {
        (*corr)[to].push_back((iter->second)[i]);
      }
      (iter->second).clear();
      (*corr).erase(iter);
    }

    // also track the poles
    if (is_medially_centered)
    {
      Point pole0 = Point(to_double((*cell_dual)[(*poles)[id0]].x()),
                          to_double((*cell_dual)[(*poles)[id0]].y()),
                          to_double((*cell_dual)[(*poles)[id0]].z()));
      Point pole1 = Point(to_double((*cell_dual)[(*poles)[id1]].x()),
                          to_double((*cell_dual)[(*poles)[id1]].y()),
                          to_double((*cell_dual)[(*poles)[id1]].z()));
      Point p1 = boost::get(*hg_point_pmap, v1);
      double dis_to_pole0 = std::sqrt(squared_distance(pole0, p1));
      double dis_to_pole1 = std::sqrt(squared_distance(pole1, p1));
      if (dis_to_pole0 < dis_to_pole1)
      {
        (*poles)[id1] = (*poles)[id0];
      }
      std::map<int, int>::iterator pole_iter = (*poles).find(id0);
      (*poles).erase(pole_iter);
    }
  }

  TriangleMeshPointPMap* hg_point_pmap;

  std::map<int, std::vector<int> >* corr;
  int max_id;

  bool is_medially_centered;
  std::map<int, int>* poles;
  std::vector<Point>* cell_dual;
};

} //namespace internal
} //namespace CGAL

/// @endcond

#endif //CGAL_MCFSKEL_TRACK_CORRESPONDENCE_H
