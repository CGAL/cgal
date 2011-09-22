// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
#define TEST_NAME    "Polyhedron's Graph Traits"
#define TEST_PROGRAM "halfedge_graph_traits_Polyhedron_3"

#include "test_prefix_Polyhedron_3.cpp"

template<class Graph>
bool test_aux ( Graph& aG )
{
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor ;
  typedef typename graph_traits<Graph>::edge_descriptor   edge_descriptor ;
  typedef typename graph_traits<Graph>::vertex_iterator   vertex_iterator ;
  typedef typename graph_traits<Graph>::edge_iterator     edge_iterator ;
  typedef typename graph_traits<Graph>::in_edge_iterator  in_edge_iterator ;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator ;

  typedef typename graph_traits<Graph const>::vertex_descriptor vertex_const_descriptor ;
  typedef typename graph_traits<Graph const>::edge_descriptor   edge_const_descriptor ;
  
  typedef typename halfedge_graph_traits<Graph>::undirected_edge_iterator undirected_edge_iterator ;
  
  bool result = false ;
  
  try
  {
    // Basic facts
    CHECK_EQUAL( num_vertices(aG), aG.size_of_vertices () ) ;
    CHECK_EQUAL( num_edges   (aG), aG.size_of_halfedges() ) ;
    
    // vtable is used to check that each vertex is reached only once by vertex_iterator
    vector<bool> vtable(aG.size_of_vertices());
    
    // Check operations on vertices
    vertex_iterator vb,ve ;
    for ( boost::tie(vb,ve) = vertices(aG) ; vb != ve ; ++ vb )
    {
      vertex_descriptor v = *vb ;
      
      // Checks that 'v' has not been reached before
      CHECK_EQUAL(vtable[v->id()],false); vtable[v->id()] = true ;

      // Degree
      CHECK_EQUAL( degree    (v,aG), v->vertex_degree() * 2) ;
      CHECK_EQUAL( out_degree(v,aG), v->vertex_degree() ) ;
      CHECK_EQUAL( in_degree (v,aG), v->vertex_degree() ) ;
    
      // Checks incident edges.
      // The number of in and out edges is counted and then asserted they match.
      size_t iec = 0, oec = 0 ;
      
      in_edge_iterator ieb,iee ;
      for ( boost::tie(ieb,iee) = in_edges(v,aG) ; ieb != iee ; ++ ieb )
      {
        edge_descriptor ie = *ieb ;
      
        // The in-edge is OK if it is incident upon 'v'
        CHECK_EQUAL(ie->vertex(),v);
        
        vertex_descriptor s = source(ie,aG);
        vertex_descriptor t = target(ie,aG);
        
        // 't' is OK if it is 'v' (since this is an in-edge)
        CHECK_EQUAL(t,v);  
        
        // 's' is OK if the outgoing edge 'ie->opposite()' is incident upon it.
        CHECK_EQUAL(s, ie->opposite()->vertex());        
        
        ++ iec ;
      }
    
      out_edge_iterator oeb,oee ;
      for ( boost::tie(oeb,oee) = out_edges(v,aG) ; oeb != oee ; ++ oeb )
      {
        edge_descriptor oe = *oeb ;
      
        // The out-edge is OK if it's opposite is incident upon 'v'
        CHECK_EQUAL(oe->opposite()->vertex(),v);
        
        vertex_descriptor s = source(oe,aG);
        vertex_descriptor t = target(oe,aG);

        // 's' is OK if it is 'v' (since this is an out-edge)
        CHECK_EQUAL(s,v);        
        
        // 't' is OK if the outgoing edge 'ie->opposite()' is incident upon it.
        CHECK_EQUAL(t, oe->vertex());        
        
        ++ oec ;
      }
      
      // Checks that the number if incoming and outgoing edges matches
      CHECK_EQUAL(iec,oec);
      
      // Checks that the number of in/out edges reached by the iterator is correct
      CHECK_EQUAL(iec,v->vertex_degree());
    }

    // Now check that all the vertices in aG are contained in the range returned by the call to vertices()
    for ( typename Graph::Vertex_const_iterator vit = aG.vertices_begin() ; vit != aG.vertices_end() ; ++ vit )
      CHECK(vtable[vit->id()]);
      
    // 'etable' is used to check that each halfedge is reached only once by edge_iterator
    vector<bool> etable(aG.size_of_halfedges());
    
    // Check operations on halfedges
    size_t ec = 0 ;
    edge_iterator eb,ee ;
    for ( boost::tie(eb,ee) = edges(aG) ; eb != ee ; ++ eb )
    {
      edge_descriptor e    = *eb ;
      edge_descriptor oe   = opposite_edge(e,aG);
      edge_descriptor ne   = next_edge(e,aG);
      edge_descriptor pe   = prev_edge(e,aG);
      edge_descriptor ccwe = next_edge_ccw(e,aG);
      edge_descriptor cwe  = next_edge_cw(e,aG);
      
      // Checks that 'e' has not been reached before
      CHECK_EQUAL(etable[e->id()],false); etable[e->id()] = true ;

      // Checks neighbors are OK      
      CHECK_EQUAL(oe,e->opposite());
      CHECK_EQUAL(ne,e->next());
      CHECK_EQUAL(pe,e->prev());
      CHECK_EQUAL(ccwe,e->opposite()->prev());
      CHECK_EQUAL(cwe,e->next()->opposite());
      
      vertex_descriptor s  = source(e,aG);
      vertex_descriptor t  = target(e,aG);
      vertex_descriptor os = source(oe,aG);
      vertex_descriptor ot = target(oe,aG);
 
      // Checks that 'o' and 'oe' are in fact opposite
      CHECK_EQUAL(s,ot);        
      CHECK_EQUAL(t,os);        
      
      // 's' is OK if the opposite edge 'oe' is incident upon it.
      CHECK_EQUAL(s, oe->vertex());        
      
      // 't' is OK if the the edge 'e' is incident upon it.
      CHECK_EQUAL(t, e->vertex());        
      
      ++ ec ;
    }
      
    // Checks that the number of halfedges reached is correct
    CHECK_EQUAL(ec,aG.size_of_halfedges());
    
    // Now check that all the halfedges in aG are contained in the range returned by the call to edges()
    for ( typename Graph::Halfedge_const_iterator eit = aG.halfedges_begin() ; eit != aG.halfedges_end() ; ++ eit )
      CHECK(etable[eit->id()]);

    // 'uetable' is used to check that only one halfedge out of each pair is reached by undirected_edge_iterator
    vector<bool> uetable(aG.size_of_halfedges());
    
    size_t uec = 0 ;
    undirected_edge_iterator ueb,uee ;
    for ( boost::tie(ueb,uee) = undirected_edges(aG) ; ueb != uee ; ++ ueb )
    {
      edge_descriptor ue = *ueb ;
      edge_descriptor oe = opposite_edge(ue,aG);
      
      // Checks that none of 'ue,oe' has been reached before
      // This also checks that 'ueb' iterates over just one halfedge in the opposing pair.
      CHECK_EQUAL(uetable[ue->id()],false); uetable[ue->id()] = true ;
      CHECK_EQUAL(uetable[oe->id()],false); uetable[oe->id()] = true ;
      
      ++ uec ;
    }

    // Checks that the number of undirected edges is exactly half the number of directed edges    
    CHECK_EQUAL(uec*2,ec);
    
    result = true ;
  }
  catch(...) {}
    
  return result ;  
}

bool test ( Polyhedron& aG )
{
  bool r = true ;
  
  r = r && test_aux<Polyhedron const>(aG);
  r = r && test_aux<Polyhedron      >(aG);
  
  return r ;
}


#include "test_suffix_Polyhedron_3.cpp"

int main( int argc, char const* argv[] )
{
  return aux_main(argc, argv);
}
