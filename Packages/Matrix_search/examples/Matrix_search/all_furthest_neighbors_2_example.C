// Copyright (c) 1998-2003  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of an example program for CGAL.  This example
// program may be used, distributed and modified without limitation.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/all_furthest_neighbors_2.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <iostream>
#include <vector>

typedef double                                    FT;

struct Kernel : public CGAL::Cartesian<FT> {};

typedef Kernel::Point_2                           Point;
typedef std::vector<int>                          Index_cont;
typedef CGAL::Polygon_2<Kernel>                   Polygon;
typedef CGAL::Random_points_in_square_2<Point>    Generator;
typedef CGAL::Ostream_iterator<int,std::ostream>  Oiterator;

int main()
{
  // generate random convex polygon:
  Polygon p;
  CGAL::random_convex_set_2(10, std::back_inserter(p), Generator(1));

  // compute all furthest neighbors:
  CGAL::all_furthest_neighbors_2(p.vertices_begin(), p.vertices_end(),
                                 Oiterator(std::cout));
  std::cout << std::endl;

  return 0;
} 
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

