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

#define CGAL_CHECK_EXPENSIVE

#include <CGAL/Real_timer.h>
#include <CGAL/Simple_cartesian.h>

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

typedef Polyhedron_3<Kernel,Polyhedron_items_with_id_3> Surface; 

typedef Surface::Vertex                                  Vertex;
typedef Surface::Vertex_iterator                         Vertex_iterator;
typedef Surface::Vertex_handle                           Vertex_handle;
typedef Surface::Vertex_const_handle                     Vertex_const_handle;
typedef Surface::Halfedge_handle                         Halfedge_handle;
typedef Surface::Halfedge_const_handle                   Halfedge_const_handle;
typedef Surface::Edge_iterator                           Edge_iterator;
typedef Surface::Facet_iterator                          Facet_iterator;
typedef Surface::Halfedge_around_vertex_const_circulator HV_circulator;
typedef Surface::Halfedge_around_facet_circulator        HF_circulator;
typedef Surface::size_type                               size_type ;


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

namespace SMS = CGAL::Surface_mesh_simplification ;

typedef Vertex_is_fixed_property_map_always_false<Surface> Vertex_is_fixed_map ;

enum Method { LT, MP } ;
enum Cache  { None, Cost, CostAndPlacement } ;

char const* method_to_string( Method aMethod )
{
  switch(aMethod)
  {
    case LT: return "LT" ; break ;
    case MP: return "MP" ; break ;
  }
  
  return "<unknown>" ;
}

char const* cache_to_string( Cache aCache )
{
  switch(aCache)
  {
    case None             : return "NoCache" ; break ;
    case Cost             : return "CostCache" ; break ;
    case CostAndPlacement : return "CostAndPlacementCache" ; break ;
  }
  
  return "<unknown>" ;
}

#include VISITOR
#include STRATEGY

typedef SMS::LindstromTurk_params LT_params ;
typedef char                      Dummy_params ;

typedef SMS::Cost_cache<Surface>                Cost_cache ;
typedef SMS::Cost_and_placement_cache<Surface>  Cost_placement_cache ;

typedef SMS::Cached_cost       <Surface>  Cached_cost ;
typedef SMS::Edge_length_cost  <Surface>  MP_cost ;
typedef SMS::LindstromTurk_cost<Surface>  LT_cost ; 
         
typedef SMS::Cached_placement<Surface>        Cached_placement ;
typedef SMS::Midpoint_placement<Surface>      MP_placement ;
typedef SMS::LindstromTurk_placement<Surface> LT_placement ;

typedef SMS::Set_no_cache<Surface> Set_no_cache ;

typedef SMS::Set_cost_cache<Surface,MP_cost>       Set_cost_cache_MP ;
typedef SMS::LindstromTurk_set_cost_cache<Surface> Set_cost_cache_LT ;

typedef SMS::Set_cost_and_placement_cache<Surface,MP_cost,MP_placement> Set_cost_placement_cache_MP ;
typedef SMS::LindstromTurk_set_cost_and_placement_cache<Surface>        Set_cost_placement_cache_LT ;

bool Test ( string aName, Method aMethod, Cache aCache )
{
  bool rSucceeded = false ;
  
  string off_name = aName ;
  
  string audit_name = aName.substr(0,aName.find_last_of(".")) + "_" + method_to_string(aMethod) + ".audit" ;
  
  ifstream off_is(off_name.c_str());
  if ( off_is )
  {
    Surface lSurface; 
    off_is >> lSurface ;
    if ( lSurface.is_valid() )
    {
      if ( lSurface.is_pure_triangle() )
      {
        ofstream audit_s(audit_name.c_str());
        if ( audit_s )
        {
          set_halfedgeds_items_id(lSurface);
          
          SMS::Count_stop_predicate<Surface> stop(1);
          
          Cached_cost get_cached_cost ;          
          MP_cost     get_MP_cost;
          LT_cost     get_LT_cost;
          
          Cached_placement get_cached_placement ;          
          MP_placement     get_MP_placement;
          LT_placement     get_LT_placement;
                    
          Set_no_cache set_no_cache ;
          
          Set_cost_cache_MP set_cost_cache_MP(get_MP_cost) ;
          Set_cost_cache_LT set_cost_cache_LT ;
          
          Set_cost_placement_cache_MP set_cost_placement_cache_MP(get_MP_cost,get_MP_placement) ;
          Set_cost_placement_cache_LT set_cost_placement_cache_LT ;
          
          Visitor lVisitor(audit_s) ;
      
          int r = -1 ;
          
          Real_timer t ; t.start();    
          switch( aMethod )
          {
            case MP:  
            
              switch ( aCache )
              {
                case None :
                
                  r = edge_collapse(lSurface
                                   ,stop
                                   ,set_cache(set_no_cache)
                                   .get_cost(get_MP_cost)
                                   .get_placement(get_MP_placement)
                                   .visitor(&lVisitor)
                                   );
                  break ;                 
                                   
                case Cost :
                
                  r = edge_collapse(lSurface
                                   ,stop
                                   ,set_cache(set_cost_cache_MP)
                                   .get_cost(get_cached_cost)
                                   .get_placement(get_MP_placement)
                                   .visitor(&lVisitor)
                                   );
                  break ;                  
                  
                case CostAndPlacement :
  
                  r = edge_collapse(lSurface
                                   ,stop
                                   ,set_cache(set_cost_placement_cache_MP)
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
                
                  r = edge_collapse(lSurface
                                   ,stop
                                   ,set_cache(set_no_cache)
                                   .get_cost(get_LT_cost)
                                   .get_placement(get_LT_placement)
                                   .visitor(&lVisitor)
                                   );
                  break ;                 
                                   
                case Cost :
                
                  r = edge_collapse(lSurface
                                   ,stop
                                   ,set_cache(set_cost_cache_LT)
                                   .get_cost(get_cached_cost)
                                   .get_placement(get_LT_placement)
                                   .visitor(&lVisitor)
                                   );
                  break ;                  
                  
                case CostAndPlacement :
                
                  r = edge_collapse(lSurface
                                   ,stop
                                   ,set_cache(set_cost_placement_cache_LT)
                                   .get_cost(get_cached_cost)
                                   .get_placement(get_cached_placement)
                                   .visitor(&lVisitor)
                                   );
                  break ;                  
                                   
              }
              
              break ;
          }
          t.stop();
          rSucceeded = true ;
        }
        else cerr << "Unable to open audit file: " << audit_name << endl ;
      }
      else
      {
        cerr << "Surface is not triangulated (has faces with more than 3 sides): " << aName << endl ;
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
  
  return rSucceeded ;
}

int main( int argc, char** argv ) 
{
  set_error_handler  (error_handler);
  set_warning_handler(error_handler);
  
  vector<string> lCases ;
  
  for ( int i = 1 ; i < argc ; ++i )
  {
    string c(argv[i]);
    if ( c.find(".off") != string::npos )
      lCases.push_back(c);
  }
   
  if ( lCases.size() == 0 )
  {
    cout << "collapse_edge_test file0 file1 ... fileN\n" ;
    return 1 ;
  } 
  else
  {
    unsigned lOK = 0 ;
    for ( vector<string>::const_iterator it = lCases.begin(); it != lCases.end() ; ++ it )
    {
     if ( Test(*it, sMethod, sCache) )
       ++ lOK ;
    }  
      
    cout << endl
         << lOK                    << " cases succedded." << endl
         << (lCases.size() - lOK ) << " cases failed." << endl ;
            
    return lOK == lCases.size() ? 0 : 1 ;
  }
}

// EOF //
