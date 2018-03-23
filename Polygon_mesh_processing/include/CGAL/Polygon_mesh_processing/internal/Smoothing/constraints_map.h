// Copyright (c) 2018 GeometryFactory (France).
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
// Author(s)     : Konstantinos Katrioplas (konst.katrioplas@gmail.com)

#ifndef CGAL_POLYGON_MESH_PROCESSING_CONSTRAINTS_MAP_H
#define CGAL_POLYGON_MESH_PROCESSING_CONSTRAINTS_MAP_H

#include <CGAL/property_map.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template<typename Descriptor>
struct Constrained_vertices_map
{
  typedef Descriptor                          key_type;
  typedef bool                                value_type;
  typedef value_type&                         reference;
  typedef boost::read_write_property_map_tag  category;

  boost::shared_ptr<std::set<Descriptor> > constrained_data;

public:
  Constrained_vertices_map() : constrained_data(new std::set<Descriptor>) {}

  friend bool get(const Constrained_vertices_map& map, const key_type& d)
  {
    typename std::set<Descriptor>::iterator it = map.constrained_data->find(d);
    return it != map.constrained_data->end() ? true : false;
  }

  friend void put(Constrained_vertices_map& map, const key_type& d)
  {
    map.constrained_data->insert(d);
  }
};

template<typename PolygonMesh, typename VertexPointMap,
         typename CotangentValue = CGAL::internal::Cotangent_value_Meyer<PolygonMesh, VertexPointMap> >
struct Edge_cotangent_weight : CotangentValue
{
    Edge_cotangent_weight(PolygonMesh& pmesh_, VertexPointMap vpmap_)
      : CotangentValue(pmesh_, vpmap_)
    {}

    PolygonMesh& pmesh()
    {
      return CotangentValue::pmesh();
    }

    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor     vertex_descriptor;

    double operator()(halfedge_descriptor he)
    {
      if(is_border_edge(he, pmesh()))
      {
        halfedge_descriptor h1 = next(he, pmesh());
        vertex_descriptor vs = source(he, pmesh());
        vertex_descriptor vt = target(he, pmesh());
        vertex_descriptor v1 = target(h1, pmesh());
        return (CotangentValue::operator ()(vs, v1, vt));
      }
      else
      {
        halfedge_descriptor h1 = next(he, pmesh());
        halfedge_descriptor h2 = prev(opposite(he, pmesh()), pmesh());
        vertex_descriptor vs = source(he, pmesh());
        vertex_descriptor vt = target(he, pmesh());
        vertex_descriptor v1 = target(h1, pmesh());
        vertex_descriptor v2 = source(h2, pmesh());
        return ( CotangentValue::operator()(vs, v1, vt) + CotangentValue::operator()(vs, v2, vt) ) / 2.0; }
    }
};

}
}
}

#endif //CGAL_POLYGON_MESH_PROCESSING_CONSTRAINTS_MAP_H
