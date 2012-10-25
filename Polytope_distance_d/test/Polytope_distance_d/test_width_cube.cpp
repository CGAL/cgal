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
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Thomas Herrmann, Lutz Kettner
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
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Width_default_traits_3.h>
#include <CGAL/Width_3.h>
#include <iostream>
#include <vector>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_real.h>
typedef leda_real                             RT;
#else
typedef double                                RT;
#endif

typedef CGAL::Cartesian<RT>                   Kernel;
typedef Kernel::RT                            RT;
typedef Kernel::Point_3                       Point_3;
typedef Kernel::Plane_3                       Plane_3;
typedef CGAL::Width_default_traits_3<Kernel>  Width_traits;
typedef CGAL::Width_3<Width_traits>           Width;
typedef CGAL::Polyhedron_3<Kernel>            Polyhedron;

int main() {
    // Create a cube using exact Cartesian coordinates
  std::vector<Point_3> points;
  points.push_back( Point_3( 1, 1, 4));
  points.push_back( Point_3(-1, 1, 4));
  points.push_back( Point_3( 1,-1, 4));
  points.push_back( Point_3(-1,-1, 4));
  points.push_back( Point_3( 1, 1, 2));
  points.push_back( Point_3(-1, 1, 2));
  points.push_back( Point_3( 1,-1, 2));
  points.push_back( Point_3(-1,-1, 2));

  // Compute convex hull of points (cube)
  Polyhedron P;
  CGAL::convex_hull_3( points.begin(), points.end(), P);

  // Compute width of cube
  Width cube(P);

  // Output of square width, width-planes, and optimal direction
  RT wnum, wdenom;
  cube.get_squared_width( wnum, wdenom);
  std::cout << "Squared Width: " << wnum << "/" << wdenom << std::endl;

  std::cout << "Direction: " << cube.get_build_direction() << std::endl;

  Plane_3 e1, e2;
  cube.get_width_planes(e1,e2);
  std::cout << "Planes: E1: " << e1 << ".  E2: " << e2 << std::endl;

  int nos = cube.get_number_of_optimal_solutions();
  std::cout << "Number of optimal solutions: "<< nos << std::endl;
  return(0);
}





