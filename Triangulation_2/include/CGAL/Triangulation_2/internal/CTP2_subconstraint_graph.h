// Copyright (c) 2020  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CTP2_SUBCONSTRAINT_GRAPH_H
#define CGAL_CTP2_SUBCONSTRAINT_GRAPH_H

#include <CGAL/license/Triangulation_2.h>

#include <CGAL/boost/graph/internal/graph_traits_2D_triangulation_helper.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>

#include <CGAL/property_map.h>

namespace CGAL
{

namespace internal
{

template <typename CTP2>
class CTP2_subconstraint_graph
{
  CTP2& ctp2;
public:

  typedef typename CTP2::Vertex_handle vertex_descriptor;
  typedef typename CTP2::Subconstraint edge_descriptor;
  typedef boost::undirected_tag directed_category;
  typedef boost::disallow_parallel_edge_tag edge_parallel_category;
  struct CTP2_graph_traversal_category :
    public virtual boost::bidirectional_graph_tag,
    public virtual boost::adjacency_graph_tag,
    public virtual boost::edge_list_graph_tag,
    public virtual boost::vertex_list_graph_tag
  { };
  typedef CTP2_graph_traversal_category traversal_category;
  typedef internal::Dereference_to_handle_enforcer<
            CTP2,
            typename CTP2::Finite_vertices_iterator,
            vertex_descriptor>                                                vertex_iterator;

  typedef typename CTP2::Subconstraint_iterator::value_type Subconstr_it_v_t;
  typedef First_of_pair_property_map<Subconstr_it_v_t> Subconstr_map;
  typedef Property_map_to_unary_function<Subconstr_map> Subconstr_uf;
  typedef boost::transform_iterator<Subconstr_uf, typename CTP2::Subconstraint_iterator> edge_iterator;

  CTP2_subconstraint_graph (CTP2& ctp2) : ctp2(ctp2) { }

  friend Iterator_range<vertex_iterator> vertices (const CTP2_subconstraint_graph& g)
  {
    return make_range (vertex_iterator(g.ctp2.finite_vertices_begin()),
                       vertex_iterator(g.ctp2.finite_vertices_end()));
  }

  friend Iterator_range<edge_iterator> edges (const CTP2_subconstraint_graph& g)
  {
    return make_range (boost::make_transform_iterator(g.ctp2.subconstraints_begin(), Subconstr_uf(Subconstr_map())),
                       boost::make_transform_iterator(g.ctp2.subconstraints_end(), Subconstr_uf(Subconstr_map())));
  }

  friend vertex_descriptor source (edge_descriptor ed, const CTP2_subconstraint_graph&)
  {
    return ed.first;
  }

  friend vertex_descriptor target (edge_descriptor ed, const CTP2_subconstraint_graph&)
  {
    return ed.second;
  }
};


template <typename CTP2>
class CTP2_graph_visitor
{
private:
  CTP2& ctp2;
  std::vector<typename CTP2::Constraint_id> to_remove;
  typename CTP2::Constraint_id current;
  typename CTP2::Vertex_handle latest_vertex;

public:

  CTP2_graph_visitor (CTP2& ctp2) : ctp2 (ctp2) { }

  void start_new_polyline()
  {
    latest_vertex = typename CTP2::Vertex_handle();
    current = typename CTP2::Constraint_id();
  }

  void add_node (typename CTP2::Vertex_handle vh)
  {
    if (latest_vertex != typename CTP2::Vertex_handle())
    {
      to_remove.push_back (ctp2.context(latest_vertex, vh).id());
      typename CTP2::Constraint_id cid = ctp2.insert_constraint(latest_vertex, vh);
      if (current == typename CTP2::Constraint_id())
        current = cid;
      else
        current = ctp2.concatenate (current, cid);
    }
    latest_vertex = vh;
  }

  void end_polyline()
  {
    for (typename CTP2::Constraint_id id : to_remove)
      ctp2.remove_constraint(id);
    to_remove.clear();
  }
};

} // namespace internal

} // namespace CGAL

#endif // CGAL_CTP2_SUBCONSTRAINT_GRAPH_H
