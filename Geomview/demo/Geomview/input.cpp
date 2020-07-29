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

// How to enter a point in Geomview.

#include <CGAL/Cartesian.h>
#include <iostream>

#ifndef CGAL_USE_GEOMVIEW
int main() {
  std::cout << "Geomview doesn't work on Windows, so..." << std::endl;
  return 0;
}
#else

#include <CGAL/IO/Geomview_stream.h>

typedef CGAL::Cartesian<double> K;

int main()
{
  CGAL::Geomview_stream gv(CGAL::Bbox_3(0, 0, 0, 350, 350, 350));

  std::cout << "Please enter a point, by right-clicking on the pickplane"
            << std::endl;

  K::Point_3 p;
  gv >> p;

  std::cout << "Here are the coordinates of the selected point : "
            << p << std::endl;

  return 0;
}
#endif
