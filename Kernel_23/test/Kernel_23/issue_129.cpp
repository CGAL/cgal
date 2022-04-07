// Copyright (c) 2015 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org)
//
// SPDX-License-Identifier: LGPL-2.1-only
//
// Author(s)     : Philipp Moeller

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <cstdlib>

// Make sure squared_radius can be called through ADL without
// triggering an access checking hard error when one of the
// internal::squared_radius functions is instantiated. This issue only
// shows up on clang with std=c++03 (the default).See
// https://github.com/CGAL/cgal/issues/129

int main()
{
  try {
    CGAL::Point_3<CGAL::Epick> a(CGAL::ORIGIN), b(a), c(a), d(a);
    CGAL::squared_radius(a, b, c, d);
  } catch(...) {}
  return EXIT_SUCCESS;
}
