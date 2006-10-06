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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Straight_skeleton_2/test/Straight_skeleton_2/test_sls_builder.cpp $
// $Id: test_sls_builder.cpp 32700 2006-07-24 22:31:02Z fcacciola $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//

#include <cstdio> 
#include <fstream>
#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Unique_hash_map.h>

using namespace boost ;
using namespace std ;
using namespace CGAL ;

typedef Simple_cartesian<double> Kernel ;
typedef Polyhedron_3<Kernel> Polyhedron ;

int sOK     = 0 ; 
int sFailed = 0 ;
 
#define CHECK(pred) \
        if (!(pred)) \
        { \
          cerr << "Assertion failure: " << #pred << endl \
               << "File:" << __FILE__ << endl \
               << "Line:" << __LINE__ << endl ; \
          throw 0 ; \
        }
        
#define CHECK_EQUAL(x,y)     CHECK(((x)==(y)))
#define CHECK_NOT_EQUAL(x,y) CHECK(((x)!=(y)))

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
    Unique_hash_map<vertex_const_descriptor,bool> vtable(false,aG.size_of_vertices());
    
    // Check operations on vertices
    vertex_iterator vb,ve ;
    for ( tie(vb,ve) = vertices(aG) ; vb != ve ; ++ vb )
    {
      vertex_descriptor v = *vb ;
      
      // Checks that 'v' has not been reached before
      CHECK_EQUAL(vtable[v],false); vtable[v] = true ;

      // Degree
      CHECK_EQUAL( degree    (v,aG), v->vertex_degree() * 2) ;
      CHECK_EQUAL( out_degree(v,aG), v->vertex_degree() ) ;
      CHECK_EQUAL( in_degree (v,aG), v->vertex_degree() ) ;
    
      // Checks incident edges.
      // The number of in and out edges is counted and then asserted they match.
      size_t iec = 0, oec = 0 ;
      
      in_edge_iterator ieb,iee ;
      for ( tie(ieb,iee) = in_edges(v,aG) ; ieb != iee ; ++ ieb )
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
      for ( tie(oeb,oee) = out_edges(v,aG) ; oeb != oee ; ++ oeb )
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
      CHECK(vtable[vit]);
      
    // 'etable' is used to check that each halfedge is reached only once by edge_iterator
    Unique_hash_map<edge_const_descriptor,bool> etable(false,aG.size_of_halfedges());
    
    // Check operations on halfedges
    size_t ec = 0 ;
    edge_iterator eb,ee ;
    for ( tie(eb,ee) = edges(aG) ; eb != ee ; ++ eb )
    {
      edge_descriptor e    = *eb ;
      edge_descriptor oe   = opposite_edge(e,aG);
      edge_descriptor ne   = next_edge(e,aG);
      edge_descriptor pe   = prev_edge(e,aG);
      edge_descriptor ccwe = next_edge_ccw(e,aG);
      edge_descriptor cwe  = next_edge_cw(e,aG);
      
      // Checks that 'e' has not been reached before
      CHECK_EQUAL(etable[e],false); etable[e] = true ;

      // Checks neighbors are OK      
      CHECK_EQUAL(oe,e->opposite());
      CHECK_EQUAL(ne,e->next());
      CHECK_EQUAL(pe,e->prev());
      CHECK_EQUAL(ccwe,e->prev()->opposite());
      CHECK_EQUAL(cwe,e->opposite()->next());
      
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
      CHECK(etable[eit]);

    // 'uetable' is used to check that only one halfedge out of each pair is reached by undirected_edge_iterator
    Unique_hash_map<edge_const_descriptor,bool> uetable(0,aG.size_of_halfedges());
    
    size_t uec = 0 ;
    undirected_edge_iterator ueb,uee ;
    for ( tie(ueb,uee) = undirected_edges(aG) ; ueb != uee ; ++ ueb )
    {
      edge_descriptor ue = *eb ;
      edge_descriptor oe = opposite_edge(ue,aG);
      
      // Checks that none of 'ue,oe' has been reached before
      CHECK_EQUAL(uetable[ue],false); uetable[ue] = true ;
      CHECK_EQUAL(uetable[oe],false); uetable[oe] = true ;
      
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

bool test( string off_file )
{
  bool rContinue = true ;
  
  std::size_t extpos = off_file.find_last_of('.') ;
  if ( extpos != string::npos && off_file.substr(extpos) == ".off" )
  {
    ifstream is(off_file.c_str());
    if ( is )
    {
      Polyhedron lPoly ;
      is >> lPoly ;
      
      bool ok = test(lPoly) ;
      
      if ( ok )
           ++ sOK ;
      else ++ sFailed ;  
      
      cout << ( ok ? "OK" : "FAILED!" ) << endl ;
    }
    else cerr << "Unable to load input .off file: " << off_file << endl ;
  }
  else cerr << "Input file must have .off extension: " << off_file << endl ;
  
  return rContinue ;
}

// This is here only to allow a breakpoint to be placed so I can trace back the problem.
void error_handler ( char const* what, char const* expr, char const* file, int line, char const* msg )
{
  cerr << "CGAL error: " << what << " violation!\n"
            << "Expr: " << expr << endl
            << "File: " << file << endl 
            << "Line: " << line << endl;
  if ( msg != 0)
      cerr << "Explanation:" << msg << endl;
      
  // Avoid an abort()    
  throw runtime_error(msg);     
}

int main( int argc, char const* argv[] )
{
  cout << "Testing Polyhedron's Graph Traits\n";
  
  CGAL::set_error_handler  (error_handler);
  CGAL::set_warning_handler(error_handler);
  
  bool print_usage = false ;
  bool nop = false ;
  
  string folder = "" ;
  
  for ( int i = 1 ; i < argc ; ++ i )  
  {
    if ( argv[i][0] == '#' )
    {
      folder = string(&argv[i][1]);
      cout << "Input folder: " << folder << endl ;
      break ;
    }  
  }
  
  vector<string> samples ;
  
  for ( int i = 1 ; i < argc ; ++ i )  
    if ( argv[i][0] != '-' && argv[i][0] != '#' )
      samples.push_back(folder+ string(argv[i])); 
  
  if ( samples.size() > 0 ) 
  {
    for ( vector<string>::const_iterator it = samples.begin() ; it != samples.end() ; ++ it )
    {
      if ( !nop )
      {
        if (!test(*it) )
          break ;
      }
      else
        cout << *it << endl ; 
    }    

    int lTotal = sOK + sFailed ;
     
    if ( lTotal > 0 )
    {
      cout << "Total cases: " << lTotal << endl
                << "Succeeded cases: " << sOK << endl
                << "Failed cases: " << sFailed << endl
                << "Failure ratio: " << ((double)sFailed/lTotal*100.0) << "%\n" ;
    } 
    
  }
  else print_usage = true ;

  if ( print_usage )
  {
    cout << "USAGE: graph_traits_Polyhedron_3 #folder file0 file1 ... fileN\n"
              << "  If file is '*' then all files with extension .off in the current folder are loaded\n" ;
  }

  return sFailed == 0 ? 0 : 1 ;
}

