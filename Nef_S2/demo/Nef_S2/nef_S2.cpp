// Copyright (c) 2004  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#include <CGAL/basic.h>

#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/Nef_S2/create_random_Nef_S2.h>
#include <CGAL/IO/Qt_widget_Nef_S2.h>
#include <qapplication.h>

typedef CGAL::Exact_integer RT;
typedef CGAL::Homogeneous<RT> Kernel;
typedef CGAL::Nef_polyhedron_S2<Kernel> Nef_polyhedron_S2;

int main(int argc, char* argv[]) {

  Nef_polyhedron_S2 S;
  create_random_Nef_S2(S,5);

  QApplication a(argc, argv);
  CGAL::Qt_widget_Nef_S2<Nef_polyhedron_S2>* w =
    new CGAL::Qt_widget_Nef_S2<Nef_polyhedron_S2>(S);
  a.setMainWidget(w);
  w->show();
  return a.exec();
}
