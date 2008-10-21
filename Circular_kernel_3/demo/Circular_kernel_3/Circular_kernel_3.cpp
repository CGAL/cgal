// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// $URL
// $Id 
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado, 
//             Julien Hazebrouck, Damien Leroy

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#include <CGAL/basic.h>

#include "Circular_kernel_3.moc"

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
