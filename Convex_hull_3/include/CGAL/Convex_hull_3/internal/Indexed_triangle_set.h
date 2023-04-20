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
#include <CGAL/Container_helper.h>
#include <boost/graph/graph_traits.hpp>

#include <vector>
#include <iostream>
#include <array>
#include <list>

namespace CGAL {
namespace Convex_hull_3 {
namespace internal {

template <typename V, typename F>
struct Indexed_triangle_set
{
  V& vertices;
  F& faces;

  typedef typename std::iterator_traits<typename F::iterator>::value_type Index_triple;
  typedef typename std::iterator_traits<typename Index_triple::iterator>::value_type Index;

  Indexed_triangle_set(V& vertices,
                       F& faces)
    : vertices(vertices), faces(faces)
  {}
};

template <typename P, typename V, typename F>
void add_isolated_points(const P& point, Indexed_triangle_set<V,F>& its)
{
  its.vertices.push_back(point);
}


template <typename P, typename V, typename F>
void copy_ch2_to_face_graph(const std::list<P>& CH_2,
                            Indexed_triangle_set<V,F>& its)
{
  its.vertices.reserve(CH_2.size());
  CGAL::internal::resize(its.faces, CH_2.size()-2);
  for(const P& p : CH_2){
    its.vertices.push_back(p);
  }

  typedef typename Indexed_triangle_set<V,F>::Index Index;

  for(std::size_t i = 1; i < CH_2.size()-1; ++i){
      CGAL::internal::resize(its.faces[i-1], 3);
      its.faces[i-1][0] = static_cast<Index>(i);
      its.faces[i-1][1] = static_cast<Index>(i + 1);
      its.faces[i-1][2] = static_cast<Index>(i + 2);
  }

}

} // namespace internal
} // namespace Convex_hull_3

template <typename TDS, typename V, typename F>
void copy_face_graph(const TDS& tds, Convex_hull_3::internal::Indexed_triangle_set<V,F>& its)
{
  typedef typename TDS::Vertex_iterator Vertex_iterator;
  typedef typename TDS::Face_iterator Face_iterator;
  CGAL::internal::resize(its.vertices, tds.number_of_vertices());
  CGAL::internal::resize(its.faces, tds.number_of_faces());
  typename Convex_hull_3::internal::Indexed_triangle_set<V,F>::Index i = 0;
  for(Vertex_iterator vit = tds.vertices_begin(); vit != tds.vertices_end(); ++vit){
      its.vertices[i] = vit->point();
    vit->info() = i++;
  }

  i = 0;
  for (Face_iterator fit = tds.faces_begin(); fit != tds.faces_end(); ++fit) {
    CGAL::internal::resize(its.faces[i], 3);
    its.faces[i][0] = fit->vertex(0)->info();
    its.faces[i][1] = fit->vertex(1)->info();
    its.faces[i][2] = fit->vertex(2)->info();
    ++i;
  }
}


template <typename V, typename F>
void clear(Convex_hull_3::internal::Indexed_triangle_set<V,F>& its)
{
  CGAL::internal::resize(its.vertices, 0);
  CGAL::internal::resize(its.faces, 0);
}


template <typename P, typename V, typename F>
void make_tetrahedron(const P& p0, const P&p1, const P& p2, const P& p3,
                      Convex_hull_3::internal::Indexed_triangle_set<V,F>& its)
{ CGAL::internal::resize(its.vertices, 4);
  CGAL::internal::resize(its.faces, 4);

  its.vertices[0] = p0;
  its.vertices[1] = p1;
  its.vertices[2] = p2;
  its.vertices[3] = p3;
  for(std::size_t i = 0; i < 4; ++i){
    CGAL::internal::resize(its.faces[i], 3);
  }
  its.faces[0][0] = 0;
  its.faces[0][1] = 1;
  its.faces[0][2] = 2;
  its.faces[1][0] = 1;
  its.faces[1][1] = 0;
  its.faces[1][2] = 3;
  its.faces[2][0] = 3;
  its.faces[2][1] = 0;
  its.faces[2][2] = 2;
  its.faces[3][0] = 2;
  its.faces[3][1] = 1;
  its.faces[3][2] = 3;
}

} // namespace CGAL

namespace boost {

// this partial specialization is needed so that the general overload
// for make_tetrahedron can be eliminated as halfedge_descriptor is
// used in the returned type
template <typename V, typename F>
  struct graph_traits<CGAL::Convex_hull_3::internal::Indexed_triangle_set<V,F>>
{
  typedef void* halfedge_descriptor;
};

} // namespace boost

#endif  // CGAL_INDEXED_TRIANGLE_SET_H
