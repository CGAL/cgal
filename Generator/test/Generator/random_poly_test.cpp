// ============================================================================
//
// Copyright (c) 2000 The GALIA Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : test/Generator/random_poly_test.C
// package       : $CGAL_Package: Generator 2.12 (28 Jul 1999) $
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Random Simple Polygons: Test Program
// ============================================================================

#define CGAL_DONT_SHUFFLE_IN_RANDOM_POLYGON_2

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Homogeneous.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Polygon_2.h>

#include <array>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <list>

template <typename K>
void test_fold()
{
  // (-5, -7) is on [(-9, 5); (-2, -16)]
  std::array<typename K::Point_2, 4> input {typename K::Point_2{-5,-7},{7,-85},{-9,5},{-2,-16}};

  int i = 0;
  do
  {
    std::cout << "permutation #" << i++ << std::endl;
    for(const auto& pt : input)
      std::cout << " (" << pt << ")";
    std::cout << std::endl;

    CGAL::Polygon_2<K> polygon;
    CGAL::random_polygon_2(input.size(), std::back_inserter(polygon), input.begin());

    if (! polygon.is_simple())
    {
      std::cerr << "ERROR: polygon is not simple." << std::endl;
      assert(false);
    }
  }
  while(std::next_permutation(input.begin(), input.end()));
}

template <typename K>
void test_random()
{
  typedef CGAL::Point_2< K >                                         Point_2;
  typedef std::list<Point_2>                                         Container;
  typedef CGAL::Polygon_2<K, Container>                              Polygon_2;
  typedef CGAL::Creator_uniform_2<double, Point_2>                   Creator;
  typedef CGAL::Random_points_in_square_2<Point_2, Creator>          Point_generator;

  Polygon_2 polygon;
  int n = 50;

  // create a polygon
  CGAL::random_polygon_2(n, std::back_inserter(polygon), Point_generator(0.5));

  // make sure it is simple
  if (! polygon.is_simple())
  {
     std::cerr << "ERROR: polygon is not simple." << std::endl;
     assert(false);
  }
}

int main()
{
  typedef CGAL::Simple_cartesian<double> CK;
  typedef CGAL::Homogeneous<double> HK;

  test_random<CK>();
  test_random<HK>();

  test_fold<CK>();
  test_fold<HK>();

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
