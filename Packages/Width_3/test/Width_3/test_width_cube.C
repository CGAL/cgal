// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : test/test_width_cube.C
// package       : Width_3 (1.6)
// chapter       : Geometric Optimisation
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Thomas Herrmann
// maintainer    : Thomas Herrmann <herrmann@ifor.math.ethz.ch>
// coordinator   : ETH Zuerich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// implementation: 3D Width of a Point Set
// ============================================================================

// short cuts for M$-VC++
#ifdef _MSC_VER
#  define Polyhedron_facet_base_3  Pfb3
#  define Vertex_max_base          Vmb
#  define Halfedge_max_base        Hmb
#endif

#include <CGAL/Cartesian.h>
#include <CGAL/leda_real.h>
#include <CGAL/Point_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <iostream>
#include <vector>
#include <CGAL/Width_default_traits_3.h>
#include <CGAL/Width_3.h>

typedef leda_real NrType;
typedef CGAL::Cartesian<NrType> RepClass;
typedef RepClass::RT RT;
typedef CGAL::Point_3<RepClass> Point;
typedef CGAL::Width_default_traits_3<RepClass> Widthtraits;
typedef CGAL::Width_3<Widthtraits> Width;
typedef CGAL::Halfedge_data_structure_polyhedron_default_3<RepClass> HDS;
typedef CGAL::Polyhedron_default_traits_3<RepClass> Polytraits;
typedef CGAL::Polyhedron_3<Polytraits,HDS> Polyhedron;


int main(int argc, char* argv[]) {
  //*** Cube: Cartesian < leda_real >***
  std::vector<Point> pointlist;
  Point p1(1,1,1);
  Point q1(-1,1,1);
  Point r1(1,-1,1);
  Point s1(-1,-1,1);
  Point p2(1,1,-1);
  Point q2(-1,1,-1);
  Point r2(1,-1,-1);
  Point s2(-1,-1,-1);
  pointlist.push_back(p1);
  pointlist.push_back(q1);
  pointlist.push_back(r1);
  pointlist.push_back(s1);
  pointlist.push_back(p2);
  pointlist.push_back(q2);
  pointlist.push_back(r2);
  pointlist.push_back(s2);

  // Compute convex hull of pointlist (cube)
  Polyhedron P;
  CGAL::convex_hull_3(pointlist.begin(),pointlist.end(),P);
  
  // Compute width of cube
  Width cube(P);
  
  // Output of width, width-planes and optimal direction
  RT WNum, WDenom;
  cube.get_squared_width(WNum,WDenom);
  std::cout<<"Squared Width: "<<WNum<<"/"<<WDenom<<std::endl;

  CGAL::Plane_3<RepClass> e1,e2;
  std::cout<<"Direction: "<<cube.get_build_direction()<<std::endl;

  CGAL::Vector_3<RepClass> dir;
  cube.get_width_planes(e1,e2);
  std::cout<<"Planes:"<<std::endl;
  std::cout<<"E1: "<<e1<<std::endl<<"E2: "<<e2<<std::endl;

  int nos=cube.get_number_of_optimal_solutions();
  std::cout<<"Number of optimal solutions: "<<nos<<std::endl;

  return(0);
}





