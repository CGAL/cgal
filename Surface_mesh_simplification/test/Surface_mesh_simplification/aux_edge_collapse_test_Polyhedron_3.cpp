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
    
  throw std::runtime_error(expr);  
}

namespace SMS = CGAL::Surface_mesh_simplification ;

typedef Vertex_is_fixed_property_map_always_false<Surface> Vertex_is_fixed_map ;

template<class T>
ostream&  operator << ( ostream& os, optional<T> const& o )
{
  if ( o )
       return os << *o ; 
  else return os << "<null>" ;
}

#define ERROR(msg) \
          { \
            cerr << "\nError: " << msg << endl \
                 << "  File:" << __FILE__ << endl \
                 << "  Line:" << __LINE__ << endl ; \
            throw runtime_error("test error"); \
          }
          
#define CHECK_MSG(pred,msg) if (!(pred)) ERROR(msg) 
         
#define CHECK(pred) CHECK_MSG(pred,#pred)
        
#define CHECK_EQUAL(x,y)       CHECK_MSG(((x)==(y)),"assertion " << #x << ":(" << x << ")==" << #y << ":(" << y << ") FAILED")
#define CHECK_NOT_EQUAL(x,y)   CHECK_MSG(((x)!=(y)),"assertion " << #x << ":(" << x << ")!=" << #y << ":(" << y << ") FAILED")

#include VISITOR

bool Test ( string aName )
{
  bool rSucceeded = false ;
  
  try
  {
    string off_name = aName ;
    
    string audit_name = aName.substr(0,aName.find_last_of(".")) + "_" + STRATEGY_ACRN + ".audit" ;
    
    ifstream off_is(off_name.c_str());
    if ( off_is )
    {
      Surface lSurface; 
      off_is >> lSurface ;
      if ( lSurface.is_valid() )
      {
        if ( lSurface.is_pure_triangle() )
        {
          cerr << "Processing " << aName << " (" << ( lSurface.size_of_halfedges() / 2 ) << " edges)" << endl ;
          
          Visitor lVisitor(audit_name) ;
          
          set_halfedgeds_items_id(lSurface);
          
          SMS::Count_stop_predicate<Surface> stop(1);
          
#include STRATEGY_POLICIES

          Real_timer t ; t.start();    
          edge_collapse(lSurface
                       ,stop
                       ,set_cache(cache)
                       .get_cost(cost)
                       .get_placement(placement)
                       .visitor(&lVisitor)
                       );
          t.stop();
                       
          rSucceeded = lSurface.is_valid() && lSurface.is_pure_triangle() ;
          
          cerr << "\r" << aName << ( rSucceeded ? " succeeded" : "FAILED" ) << ". Elapsed time=" << t.time() << " seconds." << endl ;
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
  }
  catch ( exception const& x ) 
  {
    cerr << "Exception caught: " << x.what() << endl ;
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
    string::size_type pos = c.find_last_of(".") ;
    if ( pos != string::npos )
    {
      string ext = c.substr(pos);
      if ( ext == ".off" )
        lCases.push_back(c);
    }
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
     if ( Test(*it) )
       ++ lOK ;
    }  
      
    cout << endl
         << lOK                    << " cases succedded." << endl
         << (lCases.size() - lOK ) << " cases failed." << endl ;
            
    return lOK == lCases.size() ? 0 : 1 ;
  }
}

// EOF //
