// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of an example program for CGAL.  This example
// program may be used, distributed and modified without limitation.
//
// file: examples/Generator/rcs_example.C
// --------------------------------------
// generator of random convex set
//
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>

#include <iostream>
#include <iterator>

typedef CGAL::Cartesian< double >   K;
typedef K::Point_2                  Point_2;
typedef CGAL::Random_points_in_square_2< 
     Point_2,
     CGAL::Creator_uniform_2< double, Point_2 > >
                                    Point_generator;
int main() {
  // create 500-gon and write it into a window:
  CGAL::random_convex_set_2(
            500, 
            std::ostream_iterator<Point_2>(std::cout, "\n"),
            Point_generator( 0.5));
  return 0;
}
