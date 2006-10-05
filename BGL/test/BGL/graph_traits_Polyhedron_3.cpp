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
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Simple_cartesian<double> Kernel ;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron ;

using namespace boost ;

int sOK     = 0 ; 
int sFailed = 0 ;

using namespace std ;
 
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
  
  bool result = false ;
  
  try
  {
    CHECK_EQUAL( num_vertices(aG), aG.size_of_vertices () ) ;
    CHECK_EQUAL( num_edges   (aG), aG.size_of_halfedges() ) ;
    
    vertex_iterator vb,ve ;
    for ( tie(vb,ve) = vertices(aG) ; vb != ve ; ++ vb )
    {
      vertex_descriptor v = *vb ;
      CHECK_EQUAL( degree    (v,aG), v->vertex_degree() * 2) ;
      CHECK_EQUAL( out_degree(v,aG), v->vertex_degree() ) ;
      CHECK_EQUAL( in_degree (v,aG), v->vertex_degree() ) ;
    
      in_edge_iterator ieb,iee ;
      for ( tie(ieb,iee) = in_edges(v,aG) ; ieb != iee ; ++ ieb )
      {
        edge_descriptor ie = *ieb ;
      
        vertex_descriptor s = source(ie,aG);
        vertex_descriptor t = target(ie,aG);

        CHECK_NOT_EQUAL(s,t);        
        CHECK_EQUAL    (t,v);        
        CHECK_EQUAL    (s, ie->opposite()->vertex());        
      }
    
      out_edge_iterator oeb,oee ;
      for ( tie(oeb,oee) = out_edges(v,aG) ; oeb != oee ; ++ oeb )
      {
        edge_descriptor oe = *oeb ;
      
        vertex_descriptor s = source(oe,aG);
        vertex_descriptor t = target(oe,aG);

        CHECK_NOT_EQUAL(s,t);        
        CHECK_EQUAL    (s,v);        
        CHECK_EQUAL    (t, oe->vertex());        
      }
    }

    edge_iterator eb,ee ;
    for ( tie(eb,ee) = edges(aG) ; eb != ee ; ++ eb )
    {
      edge_descriptor e = *eb ;
      
      vertex_descriptor s = source(e,aG);
      vertex_descriptor t = target(e,aG);

      CHECK_NOT_EQUAL(s,t);        
    }
      
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

