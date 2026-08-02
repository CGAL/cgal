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

int main() {
  typedef CGAL::Exact_predicates_exact_constructions_kernel     Kernel;
  typedef CGAL::Arr_segment_traits_2<Kernel>                    Traits_2;
  typedef CGAL::Arrangement_2<Traits_2>                         Arrangement_2;
  typedef CGAL::Watchtower_k_crossing_visibility_2<Arrangement_2>
                                        Watchtower_k_crossing_visibility_2;

  Watchtower_k_crossing_visibility_2 visibility;
  CGAL_USE(visibility);

  return 0;
}
