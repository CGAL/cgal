// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of an example program for CGAL.  This example
// program may be used, distributed and modified without limitation.
//
// file: examples/Generator/random_poly_example.C
// ----------------------------------------------
// program generting a random simple polygon 

#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Cartesian< double >                  K;
typedef K::Point_2                                 Point_2;
typedef std::list<Point_2>                         Container;
typedef CGAL::Polygon_2<K, Container>              Polygon_2;
typedef CGAL::Random_points_in_square_2< Point_2 > Point_generator;

int main() {
  Polygon_2 polygon;
  // create 50-gon and write it into a window:
  CGAL::random_polygon_2(50, std::back_inserter(polygon), 
                         Point_generator(0.5));
  std::cout << polygon;
  return 0;
}
