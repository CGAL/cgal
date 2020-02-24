// Copyright (c) 2015 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
namespace Mesh_3 {
namespace internal {

template <typename Graph, typename Point_3, typename NT,
          typename InterpolationFunctor>
struct Graph_manipulations
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;

  std::map<Point_3, vertex_descriptor> p2v;
  Graph& g;
  InterpolationFunctor interpolate;

  Graph_manipulations(Graph& g,
                      InterpolationFunctor interpolate = InterpolationFunctor())
    : g(g)
    , interpolate(interpolate)
  {}

  static vertex_descriptor null_vertex() {
    return boost::graph_traits<Graph>::null_vertex();
  }

  vertex_descriptor get_vertex(const Point_3& p, bool force_terminal) {
    typename std::map<Point_3, vertex_descriptor>::iterator
      it = p2v.find(p);
    if(it == p2v.end()){
      vertex_descriptor v0 = add_vertex(g);
      p2v[p] = v0;
      g[v0] = p;
      g[v0].force_terminal = force_terminal;
      return v0;
    } else {
      return it->second;
    }
  }

  template <typename Enriched_pixel, typename Null>
  vertex_descriptor split(const Enriched_pixel& epa,
                          const Enriched_pixel& epb,
                          const Null& null)
  {
    return split(epa.point,
                 epb.point,
                 epa.word,
                 epb.word,
                 null(epa.domain),
                 null(epb.domain),
                 epa.on_edge_of_the_cube,
                 epb.on_edge_of_the_cube);
  }

  vertex_descriptor split(const Point_3& a, const Point_3& b,
                          const NT v_a, const NT v_b,
                          bool a_is_outside, bool b_is_outside,
                          bool a_on_edge, bool b_on_edge)
  {
#ifdef CGAL_MESH_3_DEBUG_GRAPH_MANIPULATION
    std::cerr << "split(" << a << ", " << b << ", "
              << std::boolalpha << a_is_outside << ", "
              <<  std::boolalpha << b_is_outside << ")\n";
#endif // CGAL_MESH_3_DEBUG_GRAPH_MANIPULATION

    const Point_3 mid = interpolate(a, b, v_a, v_b);
    vertex_descriptor vmid = get_vertex(mid, a_on_edge && b_on_edge);
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
    std::cerr << "try_add_edge(" << v1 << " (" << g[v1].point
              << ", " << std::boolalpha << g[v1].force_terminal
              << "), " << v2 << " (" << g[v2].point
              << ", " << std::boolalpha << g[v2].force_terminal
              << "))\n";
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

} // namespace internal
} // namespace Mesh_3
} // namespace CGAL

#endif //CGAL_INTERNAL_MESH_3_INTERNAL_GRAPH_MANIPULATIONS
