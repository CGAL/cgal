// Copyright (c) 2011 GeometryFactory (France).
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

#ifndef CGAL_INTERNAL_COREFINEMENT_UTILS_H
#define CGAL_INTERNAL_COREFINEMENT_UTILS_H

#include <CGAL/license/Polygon_mesh_processing.h>


namespace CGAL{

namespace internal_IOP {

template <class Polyhedron>
struct Compare_unik_address{
  typedef typename Polyhedron::Halfedge_handle        Halfedge_handle;
  typedef typename Polyhedron::Halfedge_const_handle  Halfedge_const_handle;
  typedef typename Polyhedron::Halfedge               Halfedge;

  bool operator()(Halfedge_handle h1,Halfedge_handle h2) const {
    Halfedge* ph1=&(*h1) < &(*h1->opposite()) ? &(*h1) : &(*h1->opposite());
    Halfedge* ph2=&(*h2) < &(*h2->opposite()) ? &(*h2) : &(*h2->opposite());
    return  ph1 < ph2;
  }

  bool operator()(Halfedge_const_handle h1,Halfedge_const_handle h2) const {
    const Halfedge* ph1=&(*h1) < &(*h1->opposite()) ? &(*h1) : &(*h1->opposite());
    const Halfedge* ph2=&(*h2) < &(*h2->opposite()) ? &(*h2) : &(*h2->opposite());
    return  ph1 < ph2;
  }
};

template <class Polyhedron>
struct Compare_address{
  typedef typename Polyhedron::Halfedge_handle        Halfedge_handle;
  typedef typename Polyhedron::Halfedge_const_handle  Halfedge_const_handle;
  typedef typename Polyhedron::Halfedge               Halfedge;

  bool operator()(Halfedge_handle h1,Halfedge_handle h2) const {
    return  &(*h1) < &(*h2);
  }

  bool operator()(Halfedge_const_handle h1,Halfedge_const_handle h2) const {
    return  &(*h1) < &(*h2);
  }
};

template <class Polyhedron>
class Non_intersection_halfedge{
  typedef std::map< typename Polyhedron::Halfedge_const_handle,
                    std::pair<int,int>,
                    Compare_unik_address<Polyhedron>
                  >  Intersection_hedges_set;
  Intersection_hedges_set intersection_hedges_;
public:
  Non_intersection_halfedge(const Intersection_hedges_set& the_set) : intersection_hedges_(the_set){}


  bool operator()(typename Polyhedron::Halfedge_const_handle h) const
  {
    if (h->is_border_edge()) return false;
    return intersection_hedges_.find(h)==intersection_hedges_.end();
  }
};

} //end of namespace internal_IOP

namespace Corefinement{

template<class Polyhedron>
struct Dummy_edge_mark_property_map{
  typedef bool value_type;
  typedef value_type reference;
  typedef std::pair<typename Polyhedron::Halfedge_handle,Polyhedron*> key_type;
  typedef boost::read_write_property_map_tag category;

  Dummy_edge_mark_property_map(){}

  friend reference get(Dummy_edge_mark_property_map,key_type) {return false;}
  friend void put(Dummy_edge_mark_property_map,key_type,value_type) {}
};


template <class Halfedge_const_handle, class Border_halfedges_map>
int node_index_of_incident_vertex(Halfedge_const_handle h,
                                  const Border_halfedges_map& border_halfedges)
{
  //WARNING this may be expensive
  Halfedge_const_handle start=h;
  Halfedge_const_handle curr=start;
  do {
    typename Border_halfedges_map::const_iterator it_border =
      border_halfedges.find(curr );
    if (it_border!=border_halfedges.end())
      return it_border->first==curr?it_border->second.second
                                   :it_border->second.first;
    curr=curr->next()->opposite();
  } while(curr!=start);

  return -1;
}

template <class Vertex_handle, class Vertex_to_node_id>
int get_node_id(Vertex_handle vh,
                const Vertex_to_node_id& vertex_to_node_id)
{
  typename Vertex_to_node_id::const_iterator it=vertex_to_node_id.find(vh);
  if (it==vertex_to_node_id.end())
    return -1;
  return it->second;
}

} } //end of namespace CGAL::Corefinement

#endif // CGAL_INTERNAL_COREFINEMENT_UTILS_H
