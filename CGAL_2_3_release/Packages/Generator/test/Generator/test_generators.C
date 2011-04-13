// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : test_generators.C
// chapter       : $CGAL_Chapter: Geometric Object Generators $
// package       : $CGAL_Package: Generator 2.12 (28 Jul 1999) $
// source        : generators.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// Test Random Sources and Geometric Object Generators
// ============================================================================


#include <CGAL/basic.h>
#include <cstddef>
#include <cassert>
#include <vector>
#include <algorithm>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/copy_n.h>
#include <CGAL/random_selection.h>

using namespace CGAL;

typedef Cartesian<double>                R;


void test_point_generators_2() {
    typedef Point_2<R>                       Point_2;
    typedef Creator_uniform_2<double,Point_2>  Creator;

    // Create test point set.
    std::vector<Point_2> points;
    points.reserve(1000);
    Random_points_in_disc_2<Point_2,Creator>     g1( 100.0);
    CGAL::copy_n( g1, 100, std::back_inserter(points));
    Random_points_on_circle_2<Point_2,Creator>   g2( 100.0);
    Random_points_in_square_2<Point_2,Creator>   g3( 100.0);
    Random_points_on_square_2<Point_2,Creator>   g4( 100.0);
    Random_points_on_segment_2<Point_2,Creator>  g5( Point_2(-100,-100),
                                                    Point_2( 100, 100));
    Points_on_segment_2<Point_2>                g5a( Point_2( 50,-50),
                                                    Point_2(-50, 50),
                                                   50);
    CGAL::copy_n( g2, 100, std::back_inserter(points));
    CGAL::copy_n( g3, 100, std::back_inserter(points));
    CGAL::copy_n( g4, 100, std::back_inserter(points));
    CGAL::copy_n( g5,  50, std::back_inserter(points));
    CGAL::copy_n( g5a, 50, std::back_inserter(points));
    points_on_square_grid_2( 50.0, (std::size_t)100,
                             std::back_inserter(points), Creator());
    points_on_segment_2( Point_2(-100, 100), Point_2( 100,-100),
                         (std::size_t)100, std::back_inserter(points));
    random_selection( points.begin(), points.end(), (std::size_t)100,
                      std::back_inserter(points));
    random_collinear_points_2( points.begin(), points.end(),
                               (std::size_t)100,
                               std::back_inserter(points));

    // Check perturbation. Make sure that the result stays within
    // the 100 x 100 square. 10 pixel perturbation allowed.
    Random_points_in_square_2<Point_2,Creator>   g6( 90.0);
    std::vector<Point_2>::iterator i1 = points.end();
    CGAL::copy_n( g6, 100, std::back_inserter(points));
    std::vector<Point_2>::iterator i2 = points.end();
    perturb_points_2( i1, i2, 10.0);

    // Create a random permutation.
    std::random_shuffle( points.begin(), points.end(), default_random);

    CGAL_assertion( points.size() == 1000);
    for ( std::vector<Point_2>::iterator i = points.begin();
          i != points.end(); i++){
        CGAL_assertion( i->x() <=  100);
        CGAL_assertion( i->x() >= -100);
        CGAL_assertion( i->y() <=  100);
        CGAL_assertion( i->y() >= -100);
    }
}

void test_point_generators_3() {
    typedef Point_3<R>                       Point_3;
    typedef Creator_uniform_3<double,Point_3>  Creator;

    /* Create test point set. */
    std::vector<Point_3> points;
    points.reserve(400);
    Random_points_in_sphere_3<Point_3,Creator>     g1( 100.0);
    CGAL::copy_n( g1, 100, std::back_inserter(points));
    Random_points_on_sphere_3<Point_3,Creator>     g2( 100.0);
    Random_points_in_cube_3<Point_3,Creator>       g3( 100.0);
    CGAL::copy_n( g2, 100, std::back_inserter(points));
    CGAL::copy_n( g3, 100, std::back_inserter(points));
    random_selection( points.begin(), points.end(), 100,
                      std::back_inserter(points));

    CGAL_assertion( points.size() == 400);
    for ( std::vector<Point_3>::iterator i = points.begin();
          i != points.end(); i++){
        CGAL_assertion( i->x() <=  100);
        CGAL_assertion( i->x() >= -100);
        CGAL_assertion( i->y() <=  100);
        CGAL_assertion( i->y() >= -100);
        CGAL_assertion( i->z() <=  100);
        CGAL_assertion( i->z() >= -100);
    }
}

int main(){
    test_point_generators_2();
    test_point_generators_3();
    return 0;
}
// EOF //
