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
// $URL: $
// $Id: $
// 
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@gmail.com>


#include <iostream>
#include <iomanip>
#include <iostream>
#include <fstream>

#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_BGL.h>
#include <CGAL/Polyhedron_extended_BGL.h>
#include <CGAL/Polyhedron_BGL_properties.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>
#include <CGAL/IO/Polyhedron_iostream.h>

//#define CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE 4

#define CREATE_TESTCASE

void Surface_simplification_external_trace( std::string s )
{
  std::cout << s << std::endl ;
} 

int exit_code = 0 ;

#include <CGAL/Surface_mesh_simplification_vertex_pair_collapse.h>

#include <CGAL/Surface_mesh_simplification/Policies/Construct_minimal_collapse_data.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Midpoint_vertex_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/LindstromTurk.h>
#include <CGAL/Surface_mesh_simplification/Policies/Count_stop_pred.h>


template <class Refs, class Traits>
struct My_vertex : public CGAL::HalfedgeDS_vertex_base<Refs,CGAL::Tag_true,typename Traits::Point_3> 
{
  typedef CGAL::HalfedgeDS_vertex_base<Refs,CGAL::Tag_true,typename Traits::Point_3> Base ;
 
  My_vertex() {} 
  My_vertex( typename Traits::Point_3 p ) : Base(p) {}
  
  int ID; 
} ;

template <class Refs, class Traits>
struct My_halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs> 
{ 
  My_halfedge() {}
 
  int ID; 
};

template <class Refs, class Traits>
struct My_face : public CGAL::HalfedgeDS_face_base<Refs,CGAL::Tag_true,typename Traits::Plane_3>
{
  typedef CGAL::HalfedgeDS_face_base<Refs,CGAL::Tag_true,typename Traits::Plane_3> Base ;
  
  My_face() {}
  My_face( typename Traits::Plane_3 plane ) : Base(plane) {}
  
  int ID; 
};

struct My_items : public CGAL::Polyhedron_items_3 
{
    template < class Refs, class Traits>
    struct Vertex_wrapper {
        typedef My_vertex<Refs,Traits> Vertex;
    };
    template < class Refs, class Traits>
    struct Halfedge_wrapper {
        typedef My_halfedge<Refs,Traits>  Halfedge;
    };
    template < class Refs, class Traits>
    struct Face_wrapper {
        typedef My_face<Refs,Traits> Face;
    };
};

typedef CGAL::Simple_homogeneous<int>                        Kernel;
//typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Vector_3                                     Vector;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Polyhedron_3<Kernel,My_items>                  Polyhedron;

typedef Polyhedron::Vertex                                   Vertex;
typedef Polyhedron::Vertex_iterator                          Vertex_iterator;
typedef Polyhedron::Halfedge_handle                          Halfedge_handle;
typedef Polyhedron::Edge_iterator                            Edge_iterator;
typedef Polyhedron::Facet_iterator                           Facet_iterator;
typedef Polyhedron::Halfedge_around_vertex_const_circulator  HV_circulator;
typedef Polyhedron::Halfedge_around_facet_circulator         HF_circulator;

using namespace std ;

// This is here only to allow a breakpoint to be placed so I can trace back the problem.
void error_handler ( char const* what, char const* expr, char const* file, int line, char const* msg )
{
  cerr << "CGAL error: " << what << " violation!" << endl 
       << "Expr: " << expr << endl
       << "File: " << file << endl 
       << "Line: " << line << endl;
  if ( msg != 0)
    cerr << "Explanation:" << msg << endl;
}

using namespace CGAL::Triangulated_surface_mesh::Simplification ;

void Test ( char const* testcase )
{
  ifstream in(testcase);
  if ( in )
  {
    Polyhedron lP; 
    
    in >> lP ;
    
    cout << "Testing Lindstrom Turk simplification of suracewith " << (lP.size_of_halfedges()/2) << " edges..." << endl ;
  
    cout << setprecision(19) ;
  
    int lVertexID = 0 ;
    for ( Polyhedron::Vertex_iterator vi = lP.vertices_begin(); vi != lP.vertices_end() ; ++ vi )
      vi->ID = lVertexID ++ ;    
    
    int lHalfedgeID = 0 ;
    for ( Polyhedron::Halfedge_iterator hi = lP.halfedges_begin(); hi != lP.halfedges_end() ; ++ hi )
      hi->ID = lHalfedgeID++ ;    
    
    int lFacetID = 0 ;
    for ( Polyhedron::Facet_iterator fi = lP.facets_begin(); fi != lP.facets_end() ; ++ fi )
      fi->ID = lFacetID ++ ;    
  
  
    typedef LindstromTurk_collapse_data<Polyhedron> Collapse_data ;
    
    Construct_LindstromTurk_collapse_data<Collapse_data> Construct_collapse_data ;
    LindstromTurk_cost                   <Collapse_data> Get_cost ;
    LindstromTurk_vertex_placement       <Collapse_data> Get_vertex_point ;
    Count_stop_condition                 <Collapse_data> Should_stop(0);
        
    Collapse_data::Params lParams(1,1,1);
    
    int r = vertex_pair_collapse(lP,Construct_collapse_data,&lParams,Get_cost,Get_vertex_point,Should_stop);
        
    cout << "Finished...\n"
         << r << " edges removed.\n"
         << lP.size_of_vertices() << " vertices.\n"
         << (lP.size_of_halfedges()/2) << " edges.\n"
         << lP.size_of_facets() << " triangles.\n" 
         << ( lP.is_valid() ? " valid" : " INVALID!!" )
         << endl  ;
  }
  else
  {
    cerr << "Unable to open test file " << testcase << endl ;
    exit_code = 1 ;
  }              
}

void CreateTestCase()
{
  Polyhedron lP; 

  Point p( 1, 0, 0);
  Point q( 0, 1, 0);
  Point r( 0, 0, 1);
  Point s( 0, 0, 0);
  
  lP.make_tetrahedron( p, q, r, s); 
  
  ofstream out("data/sample0.off");
  if ( out )
    out << lP ;
}

int main( int argc, char** argv ) 
{
  CGAL::set_error_handler  (error_handler);
  CGAL::set_warning_handler(error_handler);
  
#ifdef CREATE_TESTCASE
  CreateTestCase();
#else
  for ( int i = 1 ; i < argc ; ++i )
    Test(argv[i]);   
#endif

  return exit_code;
}

// EOF //
