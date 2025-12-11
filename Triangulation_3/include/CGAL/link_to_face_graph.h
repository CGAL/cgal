// Copyright (c) 2017  GeometryFactory (France).
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
//


#ifndef CGAL_LINK_TO_FACE_GRAPH_H
#define CGAL_LINK_TO_FACE_GRAPH_H

#include <CGAL/license/Triangulation_3.h>


#include <boost/unordered_map.hpp>
#include <CGAL/array.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

namespace CGAL {

/*!
\ingroup PkgTriangulation3Ref

clears and fills the face graph `tm` with the <a href="https://en.wikipedia.org/wiki/Simplicial_complex#Closure.2C_star.2C_and_link">link</a> of triangulation vertex `vh`.
If `t.dimension()!=3`, nothing is done.

\tparam Triangulation_3 must be a \cgal 3D triangulation.
\tparam TriangleMesh must be a model of the concept `MutableFaceGraph`.

\param t the 3D triangulation
\param vh the vertex handle of the vertex
\param tm the triangle mesh

\param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

\cgalNamedParamsBegin
  \cgalParamNBegin{ignore_infinite_faces}
    \cgalParamDescription{if `true` and if `vh` is on the convex hull of `t`, infinite facets are ignored and a triangle mesh with a border is generated}
    \cgalParamType{Boolean}
    \cgalParamDefault{`true`}
  \cgalParamNEnd
  \cgalParamNBegin{vertex_to_vertex_output_iterator}
    \cgalParamDescription{an `OutputIterator` where pairs of input vertices and output descriptors are put upon creation in `tm`.}
    \cgalParamType{a class model of `OutputIterator` accepting
                   `std::pair<`Triangulation_3::Vertex_handle, `boost::graph_traits<TriangleMesh>::%vertex_descriptor>`}
    \cgalParamDefault{`Emptyset_iterator`}
  \cgalParamNEnd
  \cgalParamNBegin{vertex_point_map}
    \cgalParamDescription{a property map associating points to the vertices of `tm`}
    \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
                   as key type and `Triangulation_3::Point_3` as value type}
    \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
  \cgalParamNEnd
\cgalNamedParamsEnd

\returns the vertex descriptor of the triangle mesh `tm` corresponding to the infinite vertex of `t`,
         if `vh` is on the convex hull of the triangulation, and if `ignore_infinite_faces` is `false`.
         Otherwise, an arbitrary vertex descriptor of the triangle mesh `tm`.

\sa `convex_hull_3_to_face_graph()`
*/
template<class Triangulation_3,
         class TriangleMesh,
         class NamedParameters = parameters::Default_named_parameters>
typename boost::graph_traits<TriangleMesh>::vertex_descriptor
link_to_face_graph(const Triangulation_3& t,
                   typename Triangulation_3::Vertex_handle vh,
                   TriangleMesh& tm,
                   const NamedParameters& np = parameters::default_values())
{
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  vertex_descriptor nullvertex = boost::graph_traits<TriangleMesh>::null_vertex();

  if (t.dimension()!=3) return nullvertex;

  using Cell_handle = typename Triangulation_3::Cell_handle;
  using Vertex_handle = typename Triangulation_3::Vertex_handle;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  remove_all_elements(tm);
  vertex_descriptor inf;

  typedef boost::unordered_map<Vertex_handle, vertex_descriptor> Vertex_map;
  Vertex_map vertex_map;
  std::vector<Cell_handle>  cells;
  t.incident_cells(vh, std::back_inserter(cells));
  std::array<vertex_descriptor,3> face;

  typename GetVertexPointMap<TriangleMesh, NamedParameters>::type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(CGAL::vertex_point, tm));

  auto v2v_out = choose_parameter<Emptyset_iterator>(get_parameter(np, internal_np::vertex_to_vertex_output_iterator));

  const bool no_infinite_faces = choose_parameter(get_parameter(np, internal_np::ignore_infinite_faces), true);

  for(Cell_handle ch : cells){
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
          res.first->second = add_vertex(tm);
          put(vpm, res.first->second, vhj->point());
          *v2v_out++=std::make_pair(vhj, res.first->second);
          if(t.is_infinite(vhj)){
            inf = res.first->second;
          }
        }
        face[i] = res.first->second;
      }
    }
    if(!infinite_face){
      Euler::add_face(face,tm);
    }
  }
  return inf;
}

} //namespace CGAL

#endif //CGAL_LINK_TO_FACE_GRAPH_H
