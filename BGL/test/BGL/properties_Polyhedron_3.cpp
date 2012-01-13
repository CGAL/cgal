// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#define TEST_NAME    "Polyhedron's Graph Properties"
#define TEST_PROGRAM "properties_Polyhedron_3"

#include "test_prefix_Polyhedron_3.cpp"

#include "CGAL/boost/graph/properties_Polyhedron_3.h"

template<class Graph>
bool test_aux ( Graph& aG )
{
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor ;
  typedef typename graph_traits<Graph>::edge_descriptor   edge_descriptor ;
  typedef typename graph_traits<Graph>::vertex_iterator   vertex_iterator ;
  typedef typename graph_traits<Graph>::edge_iterator     edge_iterator ;

  typedef typename graph_traits<Graph const>::vertex_descriptor vertex_const_descriptor ;
  typedef typename graph_traits<Graph const>::edge_descriptor   edge_const_descriptor ;
  
  typedef typename halfedge_graph_traits<Graph>::undirected_edge_iterator undirected_edge_iterator ;
  typedef typename halfedge_graph_traits<Graph>::Point Point ;
  
  bool result = false ;
  
  try
  {
    // 
    typedef typename property_map<Graph,vertex_external_index_t>::type external_vertex_index_pmap ;
    typedef typename property_map<Graph,edge_external_index_t>::type external_edge_index_pmap ;
    
    external_vertex_index_pmap const& external_vidx_pmap = get(vertex_external_index,aG);
    external_edge_index_pmap   const& external_eidx_pmap = get(edge_external_index,aG);
    
    // Vertex propertiess
    vertex_iterator vb,ve ;
    for ( boost::tie(vb,ve) = vertices(aG) ; vb != ve ; ++ vb )
    {
      vertex_descriptor v = *vb ;

      //
      // vertex_point property      
      //
      Point& oldp = get(vertex_point,aG,v) ;
      CHECK_EQUAL(oldp,v->point());
      Point newp = midpoint(oldp,Point(ORIGIN));
      put(vertex_point,aG,v,newp);
      CHECK_EQUAL(v->point(),newp);
      // Checks that the const version also works      
      Point const& newp2 = get(vertex_point,static_cast<Graph const&>(aG),static_cast<vertex_const_descriptor>(v));
      CHECK_EQUAL(v->point(),newp2) ;
      
      //
      // vertex_index property      
      //
      CHECK_EQUAL(get(vertex_index,aG,v),v->id());
      
      //
      // vertex_external_index property      
      //
      CHECK_EQUAL(external_vidx_pmap[v],v->id());
    }

    // Edge propertiess
    edge_iterator eb,ee ;
    for ( boost::tie(eb,ee) = edges(aG) ; eb != ee ; ++ eb )
    {
      edge_descriptor e    = *eb ;
      edge_descriptor oe   = opposite_edge(e,aG);
      
      //
      // edge_weight property
      //
      double property_w = get(edge_weight,aG,e);
      double expected_w = squared_distance(e->vertex()->point(), e->opposite()->vertex()->point());
      CHECK_EQUAL(property_w,expected_w);
      
      //
      // edge_is_border property
      //
      CHECK_EQUAL(get(edge_is_border,aG,e),e->is_border());
      
      //
      // edge_index property      
      //
      CHECK_EQUAL(get(edge_index,aG,e),e->id());
      
      //
      // edge_external_index property      
      //
      CHECK_EQUAL(external_eidx_pmap[e],e->id());
    }

    result = true ;
  }
  catch(...) {}
    
  return result ;  
}

bool test ( Polyhedron& aG )
{
  return test_aux(aG);
}


#include "test_suffix_Polyhedron_3.cpp"

int main( int argc, char const* argv[] )
{
  return aux_main(argc, argv);
}
