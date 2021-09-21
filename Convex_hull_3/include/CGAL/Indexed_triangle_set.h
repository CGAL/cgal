// Copyright (c) 2021  GeometryFactory Sarl
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_INDEXED_TRIANGLE_SET_H
#define CGAL_INDEXED_TRIANGLE_SET_H

#include <CGAL/license/Convex_hull_3.h>

#include <boost/graph/graph_traits.hpp>

#include <vector>
#include <iostream>
#include <array>
#include <list>

namespace CGAL {
namespace Convex_hull_3 {
namespace internal {

template <typename P>
struct Indexed_triangle_set
{
  std::vector<P>& vertices;
  std::vector<std::array<int,3> >& faces;

  Indexed_triangle_set(std::vector<P>& vertices,
                       std::vector<std::array<int,3> >& faces)
    : vertices(vertices), faces(faces)
  {}
};




template <class P>
void add_isolated_points(const P& point, Indexed_triangle_set<P>& its)
{
  its.vertices.push_back(point);
}


template <typename P>
void copy_ch2_to_face_graph(const std::list<P>& CH_2,
                            Indexed_triangle_set<P>& its)
{
  std::cout << "copy_ch2_to_face_graph" << std::endl;
}

} // namespace internal
} // namespace Convex_hull_3

template <typename TDS, typename P>
void copy_face_graph(const TDS& tds, Convex_hull_3::internal::Indexed_triangle_set<P>& its)
{
  typedef typename TDS::Vertex_iterator Vertex_iterator;
  typedef typename TDS::Face_iterator Face_iterator;
  int i = 0;
  its.vertices.reserve(tds.number_of_vertices());
  its.faces.reserve(tds.number_of_faces());
  for(Vertex_iterator vit = tds.vertices_begin(); vit != tds.vertices_end(); ++vit){
      its.vertices.push_back(vit->point());
    vit->info() = i++;
  }
  for (Face_iterator fit = tds.faces_begin(); fit != tds.faces_end(); ++fit) {
      its.faces.push_back({fit->vertex(0)->info(), fit->vertex(1)->info(), fit->vertex(2)->info()});
  }
}


template <class P>
void clear(Convex_hull_3::internal::Indexed_triangle_set<P>& its)
{
  its.vertices.clear();
  its.faces.clear();
}

template <class P>
void make_tetrahedron(const P& p0, const P&p1, const P& p2, const P& p3,
                      Convex_hull_3::internal::Indexed_triangle_set<P>& its)
{
  CGAL_assertion(its.vertices.empty());
  its.vertices = {p0, p1, p2, p3};
  its.faces = { {0, 1, 2}, {1, 0, 3}, {3, 0, 2}, {2, 1, 3} };
}

} // namespace CGAL

namespace boost {

// this partial specialization is needed so that the general overload
// for make_tetrahedron can be eliminated as halfedge_descriptor is
// used in the returned type
template <class P>
struct graph_traits<CGAL::Convex_hull_3::internal::Indexed_triangle_set<P>>
{
  typedef void* halfedge_descriptor;
};

} // namespace boost

#endif  // CGAL_INDEXED_TRIANGLE_SET_H
