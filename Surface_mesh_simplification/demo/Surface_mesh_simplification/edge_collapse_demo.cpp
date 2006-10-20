// Copyright (c) 2002  Max Planck Institut fuer Informatik (Germany).
// All rights reserved.
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesh_simplification/test/Surface_mesh_simplification/LT_edge_collapse_test.cpp $
// $Id: LT_edge_collapse_test.cpp 32177 2006-07-03 11:55:13Z fcacciola $
// 
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@gmail.com>


#include <iostream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include <CGAL/Real_timer.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Constrained_triangulation_2.h>

#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

using namespace std ;
using namespace boost ;
using namespace CGAL ;

typedef Simple_cartesian<double> Kernel;
typedef Kernel::Vector_3         Vector;
typedef Kernel::Point_3          Point;

typedef Polyhedron_3<Kernel,Polyhedron_items_with_id_3> Polyhedron; 

typedef Polyhedron::Vertex                                  Vertex;
typedef Polyhedron::Vertex_iterator                         Vertex_iterator;
typedef Polyhedron::Vertex_handle                           Vertex_handle;
typedef Polyhedron::Vertex_const_handle                     Vertex_const_handle;
typedef Polyhedron::Halfedge_handle                         Halfedge_handle;
typedef Polyhedron::Halfedge_const_handle                   Halfedge_const_handle;
typedef Polyhedron::Edge_iterator                           Edge_iterator;
typedef Polyhedron::Facet_iterator                          Facet_iterator;
typedef Polyhedron::Halfedge_around_vertex_const_circulator HV_circulator;
typedef Polyhedron::Halfedge_around_facet_circulator        HF_circulator;
typedef Polyhedron::size_type                               size_type ;

struct Visitor
{
  Visitor ( size_t aRequested ) : mRequested(aRequested), mCollected(0) {}
  
  void OnStarted( Polyhedron& ) {} 
  
  void OnFinished ( Polyhedron& )
  { 
    cerr << "\n" << flush ;
  } 
  
  void OnStopConditionReached( Polyhedron& ) {} 
  
  void OnCollected( Halfedge_handle const& aEdge, bool aIsFixed, Polyhedron& )
  {
    ++ mCollected ;
    cerr << "\rEdges collected: " << mCollected << flush ;
  }                
  
  void OnSelected( Halfedge_handle const& aEdge, Polyhedron&, optional<double> const& aCost, size_t aInitial, size_t aCurrent )
  {
    if ( aCurrent == aInitial )
      cerr << "\n" << flush ;
        
    if ( mRequested < aInitial )
    {
      double n = aInitial - aCurrent ;  
      double d = aInitial - mRequested ; 
      cerr << "\r" << aCurrent << " " << ((int)(100.0*(n/d))) << "%" << flush ;
    }  
  }                
  
  void OnCollapsing(Halfedge_handle const& aEdge, Polyhedron&, optional<Point> const& aPlacement ) 
  {
  }                
  
  void OnNonCollapsable(Halfedge_handle const& aEdge, Polyhedron& ) 
  {
  }                
  
  size_t mRequested ;
  size_t mCollected ;
} ;

// This is here only to allow a breakpoint to be placed so I can trace back the problem.
void error_handler ( char const* what, char const* expr, char const* file, int line, char const* msg )
{
  cerr << "CGAL error: " << what << " violation!" << endl 
       << "Expr: " << expr << endl
       << "File: " << file << endl 
       << "Line: " << line << endl;
  if ( msg != 0)
    cerr << "Explanation:" << msg << endl;
    
  throw std::logic_error("");  
}

using namespace CGAL::Surface_mesh_simplification ;

char const* matched_alpha ( bool matched )
{
  return matched ? "matched" : "UNMATCHED" ; 
}

enum Method { LT, MP } ;
enum Cache  { None, Cost, CostPlacement } ;

char const* method_to_string( Method aMethod )
{
  switch(aMethod)
  {
    case LT: return "LindstromTurk" ; break ;
    case MP: return "Midpoint" ; break ;
  }
  
  return "<unknown>" ;
}

char const* cache_to_string( Cache aCache )
{
  switch(aCache)
  {
    case None          : return "None" ; break ;
    case Cost          : return "Cost" ; break ;
    case CostPlacement : return "CostPlacement" ; break ;
  }
  
  return "<unknown>" ;
}

typedef Cost_cache              <Polyhedron>  P_cost_cache ;
typedef Cost_and_placement_cache<Polyhedron>  P_cost_and_placement_cache ;

typedef Cached_cost       <Polyhedron>  P_cached_cost ;
typedef Edge_length_cost  <Polyhedron>  P_MP_cost ;
typedef LindstromTurk_cost<Polyhedron>  P_LT_cost ; 
          
typedef Cached_placement<Polyhedron>        P_cached_placement ;
typedef Midpoint_placement<Polyhedron>      P_MP_placement ;
typedef LindstromTurk_placement<Polyhedron> P_LT_placement ;

typedef Set_no_cache<Polyhedron>  P_set_no_cache ;

typedef Set_cost_cache<Polyhedron,P_MP_cost>     P_set_cost_cache_MP ;
typedef LindstromTurk_set_cost_cache<Polyhedron> P_set_cost_cache_LT ;

typedef Set_cost_and_placement_cache<Polyhedron,P_MP_cost,P_MP_placement> P_set_cost_and_placement_cache_MP ;
typedef LindstromTurk_set_cost_and_placement_cache<Polyhedron>            P_set_cost_and_placement_cache_LT ;


void Simplify ( int aStopA, int aStopR, bool aJustPrintSurfaceData, string aName, Method aMethod, Cache aCache )
{
  string off_name    = aName ;
  string result_name = aName+string(".out.off");
  
  ifstream off_is(off_name.c_str());
  if ( off_is )
  {
    Polyhedron lP; 
    
    scan_OFF(off_is,lP,true);
    
    if ( lP.is_valid() )
    {
      if ( lP.is_pure_triangle() )
      {
        if ( !aJustPrintSurfaceData )
        {
          size_t lRequestedEdgeCount ;
          if ( aStopA != -1 )
               lRequestedEdgeCount = aStopA ;
          else lRequestedEdgeCount = lP.size_of_halfedges() * aStopR / 200 ;
    
          cout << "Testing simplification of surface " << off_name 
               << " using " << method_to_string(aMethod) << " method" 
               << " with "  << cache_to_string(aCache) << " cache."  << endl ;
               
          cout << lP.size_of_facets() << " triangles." << endl 
               << (lP.size_of_halfedges()/2) << " edges." << endl 
               << lP.size_of_vertices() << " vertices." << endl 
               << (lP.is_closed() ? "Closed." : "Open." ) << endl 
               << "Requested edge count: " << lRequestedEdgeCount << endl ;
               
          cout << setprecision(19) ;
        
          set_halfedgeds_items_id(lP);
          
          P_cached_cost get_cached_cost ;          
          P_MP_cost     get_MP_cost;
          P_LT_cost     get_LT_cost;
          
          P_cached_placement get_cached_placement ;          
          P_MP_placement     get_MP_placement;
          P_LT_placement     get_LT_placement;
                    
          P_set_no_cache      set_no_cache ;
          
          P_set_cost_cache_MP set_cost_cache_MP(get_MP_cost) ;
          P_set_cost_cache_LT set_cost_cache_LT ;
          
          P_set_cost_and_placement_cache_MP set_cost_and_placement_cache_MP(get_MP_cost,get_MP_placement) ;
          P_set_cost_and_placement_cache_LT set_cost_and_placement_cache_LT ;
          
          Count_stop_predicate<Polyhedron> stop(lRequestedEdgeCount);
              
          Visitor lVisitor(lRequestedEdgeCount) ;
      
          int r = -1 ;
          
          Real_timer t ; t.start();    
          switch( aMethod )
          {
            case MP:  
            
              switch ( aCache )
              {
                case None :
                
                  r = edge_collapse(lP
                                   ,stop
                                   ,set_cache(set_no_cache)
                                   .get_cost(get_MP_cost)
                                   .get_placement(get_MP_placement)
                                   .visitor(&lVisitor)
                                   );
                  break ;                 
                                   
                case Cost :
                
                  r = edge_collapse(lP
                                   ,stop
                                   ,set_cache(set_cost_cache_MP)
                                   .get_cost(get_cached_cost)
                                   .get_placement(get_MP_placement)
                                   .visitor(&lVisitor)
                                   );
                  break ;                  
                  
                case CostPlacement :

                  r = edge_collapse(lP
                                   ,stop
                                   ,set_cache(set_cost_and_placement_cache_MP)
                                   .get_cost(get_cached_cost)
                                   .get_placement(get_cached_placement)
                                   .visitor(&lVisitor)
                                   );
                  break ;                  
                                   
              }
              
              break ;
              
            case LT: 
            
              switch ( aCache )
              {
                case None :
                
                  r = edge_collapse(lP
                                   ,stop
                                   ,set_cache(set_no_cache)
                                   .get_cost(get_LT_cost)
                                   .get_placement(get_LT_placement)
                                   .visitor(&lVisitor)
                                   );
                  break ;                 
                                   
                case Cost :
                
                  r = edge_collapse(lP
                                   ,stop
                                   ,set_cache(set_cost_cache_LT)
                                   .get_cost(get_cached_cost)
                                   .get_placement(get_LT_placement)
                                   .visitor(&lVisitor)
                                   );
                  break ;                  
                  
                case CostPlacement :
                
                  r = edge_collapse(lP
                                   ,stop
                                   ,set_cache(set_cost_and_placement_cache_LT)
                                   .get_cost(get_cached_cost)
                                   .get_placement(get_cached_placement)
                                   .visitor(&lVisitor)
                                   .edge_index_map(get(CGAL::edge_external_index,lP))
                                   );
                  break ;                  
                                   
              }
              
              break ;
          }
          t.stop();
          
          ofstream off_out(result_name.c_str(),ios::trunc);
          off_out << lP ;
          
          cout << "\nFinished...\n"
               << "Ellapsed time: " << t.time() << " seconds.\n" 
               << r << " edges removed.\n"
               << endl
               << lP.size_of_vertices() << " final vertices.\n"
               << (lP.size_of_halfedges()/2) << " final edges.\n"
               << lP.size_of_facets() << " final triangles.\n" 
               << ( lP.is_valid() ? " valid\n" : " INVALID!!\n" ) ;
        }   
        else
        {
          cout << off_name << ": " << lP.size_of_facets() << " triangles, " 
               << (lP.size_of_halfedges()/2) << " edges, "
               << lP.size_of_vertices() << " vertices, "
               << (lP.is_closed() ? "Closed." : "Open." ) << endl ;
        }  
      }
      else
      {
        cerr << "Surfaces is not triangulated (has faces with more than 3 sides): " << aName << endl ;
      }
    }
    else
    {
      cerr << "Invalid surface: " << aName << endl ;
    }
  }
  else
  {
    cerr << "Unable to open test file " << aName << endl ;
  }              
}

bool sPrintUsage = false ;

void add_case( string aCase, vector<string>& rCases )
{
  bool lAdded = false ;
  
  string::size_type pos = aCase.find_last_of(".") ;
  if ( pos != string::npos )
  {
    string ext = aCase.substr(pos);
    if ( ext == ".off" )
    {
      rCases.push_back(aCase);
      lAdded = true ;
    }
  }
    
  if ( !lAdded )
  {
    sPrintUsage = true ;  
    cerr << "Invalid input file. Only .off files are supported: " << aCase << endl ;
  }    
}

int main( int argc, char** argv ) 
{
  set_error_handler  (error_handler);
  set_warning_handler(error_handler);
  
  bool   lJustPrintSurfaceData = false ;
  int    lStopA = -1 ;
  int    lStopR = 20 ;
  Method lMethod = LT ;
  Cache  lCache  = CostPlacement ;
  string lFolder =""; 
  vector<string> lCases ;
        
  for ( int i = 1 ; i < argc ; ++i )
  {
    string opt(argv[i]);
    if ( opt[0] == '-' )
    {
      switch(opt[1])
      {
        case 'd' : lFolder = opt.substr(2); break ;
        case 'a' : lStopA = lexical_cast<int>(opt.substr(2)); break;
        case 'r' : lStopR = lexical_cast<int>(opt.substr(2)); break;
        case 'n' : lJustPrintSurfaceData = true ; break ;
        case 'm' : lMethod = (Method)lexical_cast<int>(opt.substr(2)); break ;
        case 'c' : lCache  = (Cache)lexical_cast<int>(opt.substr(2)); break ;
        
        default: 
          cerr << "Invalid option: " << opt << endl ;
          sPrintUsage = true ; 
      }
    }
  }
  
  for ( int i = 1 ; i < argc ; ++i )
  {
    string opt(argv[i]);
    if ( opt[0] == '@' )
    {
      string rspname = opt.substr(1) ;
      ifstream rsp(rspname.c_str());
      if ( rsp )
      {
        string line ;
        while ( getline(rsp,line) )
          add_case(lFolder+line,lCases);
      }
      else
      {
        cerr << "Cannot open response file: " << rspname << endl ;
        sPrintUsage = true ; 
      }
    }
    else if ( opt[0] != '-' )
    {
      add_case(lFolder+opt,lCases);
    }
  }
   
  if ( lCases.size() == 0 )
    sPrintUsage = true ;
    
  if ( sPrintUsage )
  {
    cout << "edge_collapse_demo <options> file0 file1 ... fileN @response_file" << endl 
         << "  options: " << endl
         << "    -m method                       method: 0=LindstromTurk[default] 1=Midpoint" << endl 
         << "    -c data_cache_level             level: 0=None 1=Only cost cached 2=Both cost and placement cached[default]" << endl
         << "    -d folder                       Specifies the folder where the files are located. " << endl 
         << "    -a absolute_max_edge_count      Sets the final number of edges as absolute number." << endl
         << "    -r relative_max_edge_count      Sets the final number of edges as a percentage." << endl
         << "    -n                              Do not simplify but simply report data of surfaces." << endl ;
    
    return 1 ;
  } 
  else
  {
    for ( vector<string>::const_iterator it = lCases.begin(); it != lCases.end() ; ++ it )
      Simplify( lStopA, lStopR, lJustPrintSurfaceData, *it, lMethod, lCache) ;
      
    return 0 ;
  }
}

// EOF //
