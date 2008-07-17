// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)
//
// $URL$
// $Id$
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//             Julien Hazebrouck
//             Damien Leroy

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
