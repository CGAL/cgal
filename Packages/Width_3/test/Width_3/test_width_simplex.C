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
// file          : test/test_width_simplex.C
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

// short cuts for MIPS
#if ( _COMPILER_VERSION == 730)
#  define  Homogeneous                         Hom
#  define  Width_3                             W3
#  define  Width_default_traits_3              Wdt3
#  define  Width_polyhedron_default_traits_3   Wpdt3
#  define  Width_halfedge_default_base         Whdb
#  define  Width_facet_default_base            Wfdb
#  define  Halfedge_data_structure_using_list  Hdsul
#endif

#include <CGAL/Homogeneous.h>
#include <CGAL/leda_integer.h>
#include <CGAL/leda_real.h>
#include <CGAL/Point_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Vector_3.h>
#include <iostream>
#include <vector>
#include <CGAL/Width_default_traits_3.h>
#include <CGAL/Width_3.h>

typedef leda_integer NrType;
typedef CGAL::Homogeneous<NrType> RepClass;
typedef RepClass::RT RT;
typedef CGAL::Point_3<RepClass> Point;
typedef CGAL::Width_default_traits_3<RepClass> Widthtraits;
typedef CGAL::Width_3<Widthtraits> Width;

int main(int argc, char* argv[]) {
  // *** Simplex: Homogeneous< leda_integer > ***
  std::vector<Point> pointlist;
  Point p(2,0,0,1);
  Point q(0,1,0,1);
  Point r(0,0,1,1);
  Point s(0,0,0,1);
  pointlist.push_back(p);
  pointlist.push_back(q);
  pointlist.push_back(r);
  pointlist.push_back(s);

  // Compute width of simplex
  Width simplex(pointlist.begin(),pointlist.end());

  // Output of width, width-planes and optimal direction
  RT WNum, WDenom;
  simplex.get_squared_width(WNum,WDenom);
  std::cout<<"Squared Width: "<<WNum<<"/"<<WDenom<<std::endl;

  CGAL::Plane_3<RepClass> e1,e2;
  std::cout<<"Direction: "<<simplex.get_build_direction()<<std::endl;

  CGAL::Vector_3<RepClass> dir;
  simplex.get_width_planes(e1,e2);
  std::cout<<"Planes:"<<std::endl;
  std::cout<<"E1: "<<e1<<std::endl<<"E2: "<<e2<<std::endl;

  std::cout<<"Number of optimal solutions: "
	   <<simplex.get_number_of_optimal_solutions()<<std::endl;

  return(0);
}





