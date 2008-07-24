// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado, 
//             Julien Hazebrouck, Damien Leroy

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)
// moc_source_file: demo_circle.h

#include <CGAL/basic.h>
#ifdef CGAL_USE_QT

#include "demo_circle.moc"

int main (int argc, char** argv) {
	QApplication app(argc, argv);
	MyWindow *windows = new MyWindow(1100, 1100);

	Sphere_3 laSphere(Point_3(0.,0.,0.), 1.5*1.5);
	Plane_3 plan;
	plan = Plane_3(Point_3(1., -1.,  1.), Point_3(1.,  1.,  0.), Point_3(1., -1., -1.));
	windows->add_cercle(Circle_3(laSphere, plan), 100);
	plan = Plane_3(Point_3(-1., 1.2,  1.), Point_3( 1., 1.2, -1.), Point_3( 1., 1.2,  1.));
	windows->add_cercle(Circle_3(laSphere, plan), 100);
	plan = Plane_3(Point_3(-1., -1.,  1.), Point_3( 1., -1.,  1.), Point_3( 1., -1., -1.));
	windows->add_cercle(Circle_3(laSphere, plan), 100);
	plan = Plane_3(Point_3(-1.,  1.,  0.), Point_3( 1., -1.,  1.), Point_3( 1., -1., -1.));
	windows->add_cercle(Circle_3(laSphere, plan), 100);
	plan = Plane_3(Point_3(-1., -1.,  1.), Point_3(-1., -1.,  0.), Point_3(-1.,  1., -1.));
	windows->add_cercle(Circle_3(laSphere, plan), 100);
	plan = Plane_3(Point_3( 1.,  0.,  0.), Point_3(-1.,  1.,  0.), Point_3(-1., -1.,  0.));
	windows->add_cercle(Circle_3(laSphere, plan), 100);

	app.setMainWidget(windows);
	windows->show();
	return app.exec();
}

#else // CGAL_USE_QT not defined:

#include <iostream>

int main()
{
  std::cout << "Sorry, this demo needs QT ..." << std::endl;
  return (0);
}

#endif // CGAL_USE_QT
