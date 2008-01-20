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
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Thomas Herrmann, Lutz Kettner
// maintainer    : Thomas Herrmann <herrmann@ifor.math.ethz.ch>
// coordinator   : ETH Zuerich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// implementation: 3D Width of a Point Set
// ============================================================================

#include <CGAL/Homogeneous.h>
#include <CGAL/Width_default_traits_3.h>
#include <CGAL/Width_3.h>
#include <iostream>
#include <vector>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer                          RT;
#else
typedef long                                  RT;
#endif

typedef CGAL::Homogeneous<RT>                 Kernel;
typedef Kernel::Point_3                       Point_3;
typedef Kernel::Plane_3                       Plane_3;
typedef CGAL::Width_default_traits_3<Kernel>  Width_traits;
typedef CGAL::Width_3<Width_traits>           Width;

int main() {
    // Create a simplex using homogeneous integer coordinates
    std::vector<Point_3> points;
    points.push_back( Point_3(2,0,0,1));
    points.push_back( Point_3(0,1,0,1));
    points.push_back( Point_3(0,0,1,1));
    points.push_back( Point_3(0,0,0,1));

    // Compute width of simplex
    Width simplex( points.begin(), points.end());

    // Output of squared width, width-planes, and optimal direction
    RT wnum, wdenom;
    simplex.get_squared_width( wnum, wdenom);
    std::cout << "Squared Width: " << wnum << "/" << wdenom << std::endl;

    std::cout << "Direction: " << simplex.get_build_direction() << std::endl;

    Plane_3  e1, e2;
    simplex.get_width_planes(e1,e2);
    std::cout << "Planes: E1: " << e1 << ".  E2: " << e2 <<std::endl;

    std::cout << "Number of optimal solutions: "
              << simplex.get_number_of_optimal_solutions() << std::endl;
    return(0);
}
