// Copyright (c) 2002  Max Planck Institut fuer Informatik (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner

#include <CGAL/basic.h>

#ifdef  _MSC_VER

#include <iostream>

int main() {
  std::cout << "Geomview doesn't work on Windows, so no demo." << std::endl;
  return 0;
}
#else // can have Geomeview

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>

typedef  CGAL::Cartesian<double>               Kernel;
typedef  Kernel::Point_3                       Point;
typedef  CGAL::Polyhedron_3<Kernel>            Polyhedron;

int main() {
    Point p( 1.0, 0.0, 0.0);
    Point q( 0.0, 1.0, 0.0);
    Point r( 0.0, 0.0, 1.0);
    Point s( 0.0, 0.0, 0.0);
    Polyhedron P;
    P.make_tetrahedron( p,q,r,s);
    CGAL::Geomview_stream geo;
    geo << CGAL::green() << P;

    // wait for a mouse click.
    Point click;
    geo >> click;
    return 0;
}

#endif
