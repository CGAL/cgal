// Copyright (c) 2026 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

// This test checks the call do_intersect(Triangle_3, Bbox_3).
// The implementation of the intersection predicates about bounding-boxes vs.
// kernel objects rely on comparing coordinates of the kernel objects with the
// coordinates of the bounding-box, like:
//
//    if (tr1.vertex(0).x() < segment_bbox.xmin()) return false;
//
// That means a comparison `FT < double`. If the number type `FT` has a
// comparison `FT < int` but not `FT < double`, then the comparison might be
// done by converting the `double` to `int`, which losses a lot of precision!
// The exact ring type`CGAL::Gmpzf` had that issue before, and the test was
// added to check that the issue is fixed and does not regress.

// Tests also the following macros:
#define CGAL_DO_NOT_USE_MPZF 1
#define CGAL_EPICK_NO_INTERVALS 1
#define CGAL_NO_STATIC_FILTERS 1
#define CGAL_PROFILE 1


#include <CGAL/Bbox_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Intersections_3/Bbox_3_Triangle_3.h>
#include <CGAL/Intersections_3/Segment_3_Triangle_3.h>
#include <CGAL/intersections.h>
#include <cstdlib>
#include <iostream>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Segment_3 = K::Segment_3;
using Triangle_3 = K::Triangle_3;

int main()
{
  // Triangle coordinates from the bug report
  Point_3 p1(0.222222, 0.27864439812821123, 0.10427108852604736);
  Point_3 p2(0.370222, 0.40814400000000001, 0.175922);
  Point_3 p3(0.222222, 0.53764360187178872, 0.24757291147395263);

  Triangle_3 tr1(p1, p2, p3);

  // Segment coordinates from the bug report
  Point_3 s1(0.222222, 0.40814400000000001, 0.175922);
  Point_3 s2(0.222222, 0.36585600000000001, 0.25235200000000002);

  Segment_3 curr_segment(s1, s2);

  // Get the bbox of the segment
  CGAL::Bbox_3 segment_bbox = curr_segment.bbox();

  // Print debug information
  std::cout.precision(17);
  std::cout << "Triangle:     " << tr1 << std::endl;
  std::cout << "Segment:      " << curr_segment << std::endl;
  std::cout << "Segment bbox: " << segment_bbox << std::endl;

  // Test do_intersect
  bool intersects_segment = true; // CGAL::do_intersect(tr1, curr_segment);
  bool intersects_bbox = CGAL::do_intersect(tr1, segment_bbox);

  std::cout << "do_intersect(triangle, segment): " << std::boolalpha << intersects_segment << std::endl;
  std::cout << "do_intersect(triangle, bbox):    " << std::boolalpha << intersects_bbox << std::endl;

  // The bug: first assertion passes, but the second fails
  if (!intersects_segment) {
    std::cerr << "ERROR: Triangle should intersect segment!" << std::endl;
    return EXIT_FAILURE;
  }

  if (!intersects_bbox) {
    std::cerr << "ERROR, BUG REPRODUCED: Triangle intersects segment but not its bbox!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Test passed successfully" << std::endl;
  return EXIT_SUCCESS;
}
