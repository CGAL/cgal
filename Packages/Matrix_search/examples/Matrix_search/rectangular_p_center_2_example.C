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
#include <CGAL/point_generators_2.h>
#include <CGAL/rectangular_p_center_2.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/algorithm.h>
#include <iostream>
#include <algorithm>
#include <vector>

typedef double                                      FT;

struct Kernel : public CGAL::Cartesian<FT> {};

typedef Kernel::Point_2                             Point;
typedef std::vector<Point>                          Cont;
typedef CGAL::Random_points_in_square_2<Point>      Generator;
typedef CGAL::Ostream_iterator<Point,std::ostream>  OIterator;

int main()
{
  int n = 10;
  int p = 2;
  OIterator cout_ip(std::cout);
  CGAL::set_pretty_mode(std::cout);

  Cont points;
  CGAL::copy_n(Generator(1), n, std::back_inserter(points));
  std::cout << "Generated Point Set:\n";
  std::copy(points.begin(), points.end(), cout_ip);

  FT p_radius;
  std::cout << "\n\n" << p << "-centers:\n";
  CGAL::rectangular_p_center_2(
    points.begin(), points.end(), cout_ip, p_radius, 3);
  std::cout << "\n\n" << p << "-radius = " << p_radius << std::endl;

  return 0;
} 
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

