// Copyright (c) 2014  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Sebastien Loriot
//


#ifndef CGAL_FACEGRAPH_POLYHEDRON_3_H
#define CGAL_FACEGRAPH_POLYHEDRON_3_H 1

#include <CGAL/basic.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Kernel_traits.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>

namespace CGAL {

template < class FaceGraph, class PointPMap, class HDS, bool clear_target_before = true >
class FaceGraph_to_Polyhedron_3 : public Modifier_base<HDS> {
    FaceGraph& fg;
    PointPMap ppmap;
public:
    typedef HDS  Halfedge_data_structure;
    FaceGraph_to_Polyhedron_3(const FaceGraph& src, PointPMap map)
    : fg(const_cast<FaceGraph&>(src))
    , ppmap(map)
    {}
    void operator()( HDS& target);
};

template < class FaceGraph, class PointPMap, class HDS, bool clear_target_before>
void
FaceGraph_to_Polyhedron_3<FaceGraph, PointPMap, HDS, clear_target_before>::
operator()(HDS& tgt)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;

  Cartesian_converter<
    typename Kernel_traits<typename boost::property_traits<PointPMap>::value_type>::Kernel,
    typename Kernel_traits<typename HDS::Vertex::Point>::Kernel
  > convert;

  if ( clear_target_before )
    tgt.clear();

  Polyhedron_incremental_builder_3<HDS> B(tgt);
  B.begin_surface( num_vertices(fg),
                   num_faces(fg),
                   num_halfedges(fg));
  std::map<vertex_descriptor, std::size_t> indices;
  std::size_t i=0;
  BOOST_FOREACH(vertex_descriptor vd, vertices(fg) )
  {
    B.add_vertex( convert( get(ppmap, vd) ) );
    indices[vd]=i++;
  }

  BOOST_FOREACH(face_descriptor fd, faces(fg))
  {
    B.begin_facet();
    halfedge_descriptor hd=halfedge(fd, fg), first=hd;
    do {
        B.add_vertex_to_facet( indices[target(edge(hd,fg), fg)] );
        hd=next(hd,fg);;
    } while( hd != first);
    B.end_facet();
  }
  B.end_surface();
}

} //namespace CGAL
#endif // CGAL_FACEGRAPH_POLYHEDRON_3_H //
// EOF //
