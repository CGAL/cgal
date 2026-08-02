// Copyright (c) 2026 Toronto Metropolitan Unversity (Canada).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s):  Yeganeh Bahoo Torudi <bahoo@torontomu.ca>
//             Teodor Cirstoiu <tcirstoiu@torontomu.ca>
//

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Watchtower_k_crossing_visibility_2.h>

#include <iostream>
#include <vector>

int main() {
  typedef CGAL::Exact_predicates_exact_constructions_kernel     Kernel;
  typedef Kernel::Point_2                                       Point_2;

  // x-monotone terrain, vertices in increasing x order
  std::vector<Point_2> terrain;
  terrain.emplace_back(0, 0);
  terrain.emplace_back(1, 3);
  terrain.emplace_back(2, 1);
  terrain.emplace_back(3, 4);
  terrain.emplace_back(4, 0);

  const unsigned int k = 1;

  const Point_2 w_cont = CGAL::continuous_watchtower_k_crossing_visibility_2(
                     terrain.begin(), terrain.end(), k);

  const Point_2 w_disc = CGAL::discrete_watchtower_k_crossing_visibility_2(
                     terrain.begin(), terrain.end(), k);

  std::cout << "continuous: " << w_cont << std::endl;
  std::cout << "discrete:   " << w_disc << std::endl;

  return 0;
}
