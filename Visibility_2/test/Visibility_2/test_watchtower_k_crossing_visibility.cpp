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
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Watchtower_k_crossing_visibility_2.h>

#include <iostream>
#include <vector>

int main() {
  typedef CGAL::Exact_predicates_exact_constructions_kernel     Kernel;
  typedef Kernel::Point_2                                       Point_2;
  typedef CGAL::Arr_segment_traits_2<Kernel>                    Traits_2;
  typedef CGAL::Arrangement_2<Traits_2>                         Arrangement_2;
  typedef CGAL::Watchtower_k_crossing_visibility<Arrangement_2>
                                        Watchtower_k_crossing_visibility;

  // x-monotone terrain, vertices in increasing x order
  std::vector<Point_2> terrain;
  terrain.emplace_back(0, 0);
  terrain.emplace_back(1, 3);
  terrain.emplace_back(2, 1);
  terrain.emplace_back(3, 4);
  terrain.emplace_back(4, 0);

  const unsigned int k = 1;

  const Watchtower_k_crossing_visibility watchtower(k);

  const Point_2 w_cont = watchtower.compute_continuous_watchtower(
                     terrain.begin(), terrain.end());

  const Point_2 w_disc = watchtower.compute_discrete_watchtower(
                     terrain.begin(), terrain.end());

  std::cout << "continuous: " << w_cont << std::endl;
  std::cout << "discrete:   " << w_disc << std::endl;

  return 0;
}
