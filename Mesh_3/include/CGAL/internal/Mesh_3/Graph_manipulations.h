// Copyright (c) 2015 GeometryFactory
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
// Author(s)     : Laurent Rineau

#ifndef CGAL_INTERNAL_MESH_3_INTERNAL_GRAPH_MANIPULATIONS
#define CGAL_INTERNAL_MESH_3_INTERNAL_GRAPH_MANIPULATIONS

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Kernel_traits.h>
// Assumes the point is a CGAL point.

#include <boost/graph/graph_traits.hpp>

#include <map>

namespace CGAL {
namespace internal {
namespace Mesh_3 {

template <typename Graph, typename Point_3>
struct Graph_manipulations
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;

  std::map<Point_3, vertex_descriptor> p2v;
  Graph& g;

  Graph_manipulations(Graph& g) : g(g) {}

  vertex_descriptor get_vertex(const Point_3& p) {
    typename std::map<Point_3, vertex_descriptor>::iterator
      it = p2v.find(p);
    if(it == p2v.end()){
      vertex_descriptor v0 = add_vertex(g);
      p2v[p] = v0;
      g[v0] = p;
      return v0;
    } else {
      return it->second;
    }
  }

  vertex_descriptor split(const Point_3& a, const Point_3& b,
                          bool a_is_outside, bool b_is_outside)
  {
#ifdef CGAL_MESH_3_DEBUG_GRAPH_MANIPULATION
    std::cerr << "split(" << a << ", " << b << ", "
              << std::boolalpha << a_is_outside << ", "
              <<  std::boolalpha << b_is_outside << ")\n";
#endif // CGAL_MESH_3_DEBUG_GRAPH_MANIPULATION

    typedef typename CGAL::Kernel_traits<Point_3>::Kernel K;
    typename K::Construct_midpoint_3 midpt
      = K().construct_midpoint_3_object();
    const Point_3 mid = a < b ? midpt(a, b) : midpt(b, a);
    vertex_descriptor vmid = get_vertex(mid);
    typename std::map<Point_3, vertex_descriptor>::iterator
      it_a = p2v.find(a),
      it_b = p2v.find(b);
    if(it_a != p2v.end() && it_b != p2v.end()) {
      vertex_descriptor va = it_a->second;
      vertex_descriptor vb = it_b->second;
      edge_descriptor edge;
      bool b;
      // test if the edge is already here, using add_edge
      boost::tie(edge, b) = add_edge(va, vb, g);
      remove_edge(edge, g);
      if(!b) {
        // The edge was already here.
        if(!a_is_outside) try_add_edge(va, vmid);
        if(!b_is_outside) try_add_edge(vb, vmid);
#ifdef CGAL_MESH_3_DEBUG_GRAPH_MANIPULATION
        std::cerr << " --> vmid = " << vmid << "\n";
#endif // CGAL_MESH_3_DEBUG_GRAPH_MANIPULATION
        return vmid;
      }
    }
#ifdef CGAL_MESH_3_DEBUG_GRAPH_MANIPULATION
    std::cerr << " --> vmid = " << vmid << "\n";
#endif // CGAL_MESH_3_DEBUG_GRAPH_MANIPULATION
    return vmid;
  }

  bool try_add_edge(vertex_descriptor v1, vertex_descriptor v2) {
#ifdef CGAL_MESH_3_DEBUG_GRAPH_MANIPULATION
    std::cerr << "try_add_edge(" << v1 << " (" << g[v1]
              << "), " << v2 << " (" << g[v2] << "))\n";
#endif // CGAL_MESH_3_DEBUG_GRAPH_MANIPULATION
    if(v1 != v2) {
      edge_descriptor edge;
      bool b;
      boost::tie(edge, b) = add_edge(v1, v2, g);
      return b;
    } else
      return false;
  }
}; // struct template Graph_manipulations

} // namespace Mesh_3
} // namespace internal
} // namespace CGAL

#endif //CGAL_INTERNAL_MESH_3_INTERNAL_GRAPH_MANIPULATIONS
