// Copyright (c) 1998-2003  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
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
  CGAL::random_convex_set_2(10, back_inserter(p), Generator(1));

  // compute all furthest neighbors:
  CGAL::all_furthest_neighbors_2(p.vertices_begin(), p.vertices_end(),
                                 Oiterator(std::cout));
  std::cout << std::endl;

  return 0;
} 
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

