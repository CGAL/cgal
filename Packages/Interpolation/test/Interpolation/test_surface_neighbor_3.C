// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : /test/Interpolation/test_surface_neighbors_2.C
// package       : Interpolation
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Julia Floetotto <Julia.Flototto@sophia.inria.fr>
//
// coordinator   : 
// ============================================================================
#include <CGAL/basic.h>
#include <CGAL/double.h>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// #include <CGAL/Homogeneous.h>
// #include <CGAL/MP_Float.h>

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>
#include <CGAL/aff_transformation_tags.h>

#include <CGAL/_test_surface_neighbors_3.C>


struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Delaunay_triangulation_3<K>                   Dt;

struct K2 : CGAL::Exact_predicates_exact_constructions_kernel {};
typedef CGAL::Delaunay_triangulation_3<K2>                   Dt2;

//Hierarchy with exact pred exact const. kernel:
typedef CGAL::Triangulation_vertex_base_3<K2>             Vb;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb>  Vbh;
typedef CGAL::Triangulation_data_structure_3<Vbh>        Tds;
typedef CGAL::Delaunay_triangulation_3<K2,Tds>            Dt1;
typedef CGAL::Triangulation_hierarchy_3<Dt1>              Dh;

// //Homogeneous Kernel:
// typedef CGAL::MP_Float                                  NT;
// typedef CGAL::Homogeneous<NT>                           K3;
// typedef CGAL::Delaunay_triangulation_3<K3>                   Dt3;

// Aff_transformation:
typedef CGAL::Aff_transformation_3<K2>                   Transformation;
 

int main()
{

  std::cout << "Using Exact_predicates_inexact_constructions_kernel: " << std::endl ;
  std::cout << "Testing surface_neighbors_3 on a sphere "; 
  _test_surface_neighbors_3_sphere( Dt() );
  std::cout << " done." << std::endl << std::endl;

 //  std::cout << "Using Homogeneous_kernel: " << std::endl ;
//   std::cout << "Testing surface_neighbors_3 on a sphere "; 
//   _test_surface_neighbors_3_sphere( Dt3() );
//   std::cout << " done." << std::endl << std::endl;

  std::cout << "Using Exact_predicates_exact_constructions_kernel: "
	    << std::endl;

  //AXIS ALIGNED CUBE
  Transformation identity(CGAL::IDENTITY);
  std::cout << "Testing surface_neighbors_3 on a cube with axis : " 
	    << identity(K2::Vector_3(0,0,1)) << ",  "
	    << identity(K2::Vector_3(0,1,0))<< ",  "
	    << identity(K2::Vector_3(1,0,0))
	    << std::endl;
  std::cout << " with grid sample points";
  _test_surface_neighbors_3_cube(Dt2(),identity);
  std::cout << " done." << std::endl; 
  
  std::cout << " with random sample points";
  _test_surface_neighbors_3_cube( Dt2(), identity,150, K2::FT(1e-2),false);
  std::cout << " done." << std::endl << std::endl; 
  
  //ROTATED CUBE
  Transformation  rotate(K2::RT(1),K2::RT(0),K2::RT(0),K2::RT(0),
			  K2::RT(0),K2::RT(0.9063),K2::RT(-0.42261826),K2::RT(0),
			  K2::RT(0),K2::RT(0.42261826),K2::RT(0.9063),K2::RT(0)); 
  std::cout << "Testing surface_neighbors_3 on a ROTATED cube "<< std::endl;
  std::cout << " with grid sample points";
  _test_surface_neighbors_3_cube(Dh(),rotate, 75, K2::FT(1e-2), true);
  std::cout << " done." << std::endl << std::endl; 

  
 //  //undersampled rotated cube
//   Transformation  rotate3(K2::RT(0.1),K2::RT(0.4),K2::RT(0.6),K2::RT(0),
// 			  K2::RT(0.3),K2::RT(0.5),K2::RT(0.1),K2::RT(0),
// 			  K2::RT(0.5),K2::RT(0.9),K2::RT(0.8),K2::RT(0)); 
//   std::cout << "Testing surface_neighbors_3 on an undersampled ROTATED cube "
// 	    << std::endl;
//   std::cout << " with grid sample points";
//   _test_surface_neighbors_3_cube(Dh(),rotate3,75, K2::FT(9), true);
//   std::cout << " done." << std::endl << std::endl; 

  return 0;
};

