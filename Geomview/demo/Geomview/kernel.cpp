// Copyright (c) 2000, 2001, 2004
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion

// Demo program for Geomview_stream with kernel objects.

#include <CGAL/Cartesian.h>
#include <iostream>

#ifndef CGAL_USE_GEOMVIEW
int main() {
  std::cout << "Geomview doesn't work on Windows, so..." << std::endl;
  return 0;
}
#else

#include <unistd.h>

#include <CGAL/intersections.h>
#include <CGAL/IO/Geomview_stream.h>

typedef CGAL::Cartesian<double> K;

int main()
{
  CGAL::Geomview_stream gv(CGAL::Bbox_3(0, 0, 0, 350, 350, 350));

  // gv.set_trace(true);
  gv.clear(); // remove the pickplane.

  gv << K::Point_2 (200, 100);
  gv << CGAL::blue();
  gv << K::Point_3 (200, 100, 100);
  gv << CGAL::red();
  gv << K::Segment_2 (K::Point_2(200, 100),
                      K::Point_2(300, 100));
  gv << CGAL::green();
  gv << K::Segment_3 (K::Point_3(200, 100, 100),
                      K::Point_3(300, 100, 200));
  gv << CGAL::deep_blue();
  gv << K::Sphere_3 (K::Point_3(100, 100, 100), 1000);
  gv << CGAL::violet();
  gv << K::Triangle_2 (K::Point_2(200, 200),
                       K::Point_2(220, 220),
                       K::Point_2(180, 220));
  gv << CGAL::orange();
  gv << K::Triangle_3 (K::Point_3(200, 200, 50),
                       K::Point_3(220, 220, 80),
                       K::Point_3(180, 220, 100));
  gv << CGAL::purple();
  gv << K::Tetrahedron_3 (K::Point_3(100, 100, 180),
                          K::Point_3(120,  70, 220),
                          K::Point_3(100, 100, 220),
                          K::Point_3(120, 150, 250));
  gv << CGAL::Bbox_2(10, 10, 30, 30);
  gv << CGAL::Bbox_3(10, 10, 10, 30, 30, 30);

  gv << CGAL::red();
  gv << K::Ray_2(K::Point_2(205,205), K::Point_2(500,500));
  gv << K::Ray_3(K::Point_3(250,250,250), K::Point_3(500,500,500));
  gv << K::Line_2(K::Point_2(195,195), K::Point_2(500,500));
  gv << K::Line_3(K::Point_3(150,150,150), K::Point_3(500,500,500));

  gv.look_recenter();

  std::cout << "Stopping in 1 minute" << std::endl;
  sleep(60);

  return 0;
}
#endif
