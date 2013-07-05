// Copyright (c) 1998-2003  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/distance_predicates_2.h>
#include <CGAL/all_furthest_neighbors_2.h>
#include <vector>

using std::equal;
using std::back_inserter;
using CGAL::random_convex_set_2;
using CGAL::has_smaller_distance_to_point;
using CGAL::squared_distance;
using CGAL::iterator_distance;

typedef double                                    FT;

struct Kernel : public CGAL::Cartesian<FT> {};

typedef Kernel::Point_2                           Point;
typedef std::vector<int>                          Index_cont;
typedef CGAL::Polygon_2<Kernel>                   Polygon_2;
typedef CGAL::Random_points_in_square_2<Point>    Generator;

#include <CGAL/squared_distance_2.h>
#include <CGAL/circulator.h>
#include <algorithm>

template < class RandomAccessIC, class OutputIterator >
OutputIterator
afn_brute_force(RandomAccessIC b, RandomAccessIC e, OutputIterator o)
{
  RandomAccessIC i1 = b;
  do {
    RandomAccessIC i2 = b;
    RandomAccessIC i = b;
    do {
      if (squared_distance(*i1, *i2) > squared_distance(*i1, *i))
        i = i2;
    } while (++i2 != e);
    *o++ = static_cast<int>(iterator_distance(b, i));
  } while (++i1 != e);
  return o;
} // afn_brute_force(b, e, o)

int main()
{
  int size[] = { 3, 5, 20, 101, 534 };
  for (int i = 0; i < 5; ++i) {
    int n = size[i];
    // generate random convex polygon:
    Polygon_2 p;
    CGAL::random_convex_set_2(n, std::back_inserter(p), Generator(1));
    // compute all furthest neighbors:
    Index_cont neighbors;
    CGAL::all_furthest_neighbors_2(p.vertices_begin(), p.vertices_end(),
                                   back_inserter(neighbors));
    // compute again brute force:
    Index_cont neighbors2;
    afn_brute_force(p.vertices_begin(), p.vertices_end(),
                    back_inserter(neighbors2));
    
    // and compare both results:
    assert( equal( neighbors.begin(),
                           neighbors.end(),
                           neighbors2.begin()));
  } // for (i = 0..4)

  // try also once with a random-acccess iterator:
  int n = 222;
  // generate random convex polygon:
  Polygon_2 p;
  CGAL::random_convex_set_2(n, std::back_inserter(p), Generator(1));
  Index_cont neighbors(n);
  CGAL::all_furthest_neighbors_2(p.vertices_begin(), p.vertices_end(),
                                 neighbors.begin());

  return 0;
} // int main()
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

