// Copyright (c) 2020,2025  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CTP2_SUBCONSTRAINT_GRAPH_H
#define CGAL_CTP2_SUBCONSTRAINT_GRAPH_H

#include <CGAL/license/Triangulation_2.h>

#include <CGAL/boost/graph/internal/graph_traits_2D_triangulation_helper.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>

namespace CGAL
{

namespace internal
{

template <typename CTP2>
class CTP2_subconstraint_graph
{
  CTP2& ctp2;
public:
  using vertex_descriptor = typename CTP2::Vertex_handle;
  using edge_descriptor = typename CTP2::Subconstraint;
  using directed_category = boost::undirected_tag;
  using edge_parallel_category = boost::disallow_parallel_edge_tag;
  struct CTP2_graph_traversal_category : public virtual boost::bidirectional_graph_tag,
                                         public virtual boost::adjacency_graph_tag,
                                         public virtual boost::edge_list_graph_tag,
                                         public virtual boost::vertex_list_graph_tag
  {};
  using traversal_category = CTP2_graph_traversal_category;
  using vertex_iterator =
      internal::Dereference_to_handle_enforcer<CTP2, typename CTP2::Finite_vertices_iterator, vertex_descriptor>;

  using edge_iterator = typename CTP2::Subconstraint_iterator;

  CTP2_subconstraint_graph(CTP2& ctp2)
      : ctp2(ctp2) {}

  friend Iterator_range<vertex_iterator> vertices(const CTP2_subconstraint_graph& g) {
    return make_range(vertex_iterator(g.ctp2.finite_vertices_begin()), vertex_iterator(g.ctp2.finite_vertices_end()));
  }

  friend Iterator_range<edge_iterator> edges(const CTP2_subconstraint_graph& g) {
    return g.ctp2.subconstraints();
  }

  friend vertex_descriptor source(edge_descriptor ed, const CTP2_subconstraint_graph&) { return ed.first; }

  friend vertex_descriptor target(edge_descriptor ed, const CTP2_subconstraint_graph&) { return ed.second; }
};

template <typename CTP2>
class CTP2_graph_visitor
{
private:
  using Constraint_id = typename CTP2::Constraint_id;
  using Vertex_handle = typename CTP2::Vertex_handle;
  CTP2& ctp2;
  std::vector<Constraint_id> to_remove;
  Constraint_id current{};
  Vertex_handle latest_vertex{};

public:
  CTP2_graph_visitor(CTP2& ctp2)
      : ctp2(ctp2) {}

  void start_new_polyline() {
    latest_vertex = Vertex_handle();
    current = Constraint_id();
  }

  void add_node(Vertex_handle vh) {
    if(latest_vertex != Vertex_handle()) {
      to_remove.push_back(ctp2.context(latest_vertex, vh).id());
      Constraint_id cid = ctp2.insert_constraint(latest_vertex, vh);
      if(current == Constraint_id())
        current = cid;
      else
        current = ctp2.concatenate(current, cid);
    }
    latest_vertex = vh;
  }

  void end_polyline() {
    for(Constraint_id id : to_remove)
      ctp2.remove_constraint(id);
    to_remove.clear();
  }
};

} // namespace internal

} // namespace CGAL

#endif // CGAL_CTP2_SUBCONSTRAINT_GRAPH_H
