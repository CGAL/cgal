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

#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

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

bool Test ( string aName )
{
  bool rSucceeded = false ;
  
  string off_name    = aName ;
  string audit_name  = aName+string(".audit");
  string result_name = aName+string(".out.off");
  
  ifstream off_is(off_name.c_str());
  if ( off_is )
  {
    Polyhedron lP; 
    off_is >> lP ;
    if ( lP.is_valid() )
    {
      if ( lP.is_pure_triangle() )
      {
        rSucceeded = true ;
      }
      else cerr << "ERROR: Surface is not purely triangular: " << aName << endl ;
    }  
    else cerr << "ERROR: Surface file is invalid: " << aName << endl ;
  }
  else cerr << "ERROR: Unable to open test file " << aName << endl ;
  
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
    
  if ( lCases.size() == 0)
  {
    cout << "collapse_edge_test file0 file1 ... fileN" << endl ;
    return 1 ;
  } 
  else
  {
    unsigned lOK = 0 ;
    for ( vector<string>::const_iterator it = lCases.begin(); it != lCases.end() ; ++ it )
    {
      if ( Test( *it ) )
        ++ lOK ;
    }  
      
    cout << endl
         << lOK                    << " cases succedded." << endl
         << (lCases.size() - lOK ) << " cases failed." << endl ;
            
    return lOK == lCases.size() ? 0 : 1 ;
  }
}

// EOF //
