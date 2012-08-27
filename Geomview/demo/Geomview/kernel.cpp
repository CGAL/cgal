// Copyright (c) 2000, 2001, 2004  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
  gv << CGAL::BLUE;
  gv << K::Point_3 (200, 100, 100);
  gv << CGAL::RED;
  gv << K::Segment_2 (K::Point_2(200, 100),
                      K::Point_2(300, 100));
  gv << CGAL::GREEN;
  gv << K::Segment_3 (K::Point_3(200, 100, 100),
                      K::Point_3(300, 100, 200));
  gv << CGAL::DEEPBLUE;
  gv << K::Sphere_3 (K::Point_3(100, 100, 100), 1000);
  gv << CGAL::VIOLET;
  gv << K::Triangle_2 (K::Point_2(200, 200),
                       K::Point_2(220, 220),
                       K::Point_2(180, 220));
  gv << CGAL::ORANGE;
  gv << K::Triangle_3 (K::Point_3(200, 200, 50),
                       K::Point_3(220, 220, 80),
                       K::Point_3(180, 220, 100));
  gv << CGAL::PURPLE;
  gv << K::Tetrahedron_3 (K::Point_3(100, 100, 180),
                          K::Point_3(120,  70, 220),
                          K::Point_3(100, 100, 220),
                          K::Point_3(120, 150, 250));
  gv << CGAL::Bbox_2(10, 10, 30, 30);
  gv << CGAL::Bbox_3(10, 10, 10, 30, 30, 30);

  gv << CGAL::RED;
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
