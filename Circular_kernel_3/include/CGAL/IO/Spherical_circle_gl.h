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
//             Sylvain Pion
//             Julien Hazebrouck
//             Damien Leroy

#ifndef CGAL_SPHERICAL_CIRCLE_GL__H
#define CGAL_SPHERICAL_CIRCLE_GL__H

#include <cmath>

namespace CGAL {

#define COLOR_POINT(p)            glColor3f(abs(p.x()),abs(p.y()),abs(p.z()))
#define AFFICHER_POINT_CONTOUR(p) COLOR_POINT(p),glVertex3f(CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z()))
#define TRACE_SEGMENT(p1,p2)      AFFICHER_POINT_CONTOUR(p1),AFFICHER_POINT_CONTOUR(p2)
#define COLOR_ROUGE               glColor3f(1.0, 0.0, 0.0)
#define COLOR_POINT(p)            glColor3f(abs(p.x()),abs(p.y()),abs(p.z()))
#define IMPRIMER_POINT(n,p)       std::cout << n << " = (" << p.x() << ", " << p.y() << ", " << p.z() << ")" << std::endl
#define VECTOR_LENGTH(v)          sqrt((v.x()*v.x())+(v.y()*v.y())+(v.z()*v.z()))

/**
 * Function drawing the contour of the 3D circle passed as argument with the desired precision.
 */
template <class SK>
void dessiner_spherical_circle (const typename SK::Circle_3& circle, int precision) {
	typedef typename SK::Point_3        Point_3;
	typedef typename SK::Vector_3       Vector_3;
	typedef typename SK::Has_on_3       Has_on_3;
	const Point_3& centre = circle.center();
	//IMPRIMER_POINT("centre", centre);
	Vector_3 u = circle.supporting_plane().base1();
	u = u / VECTOR_LENGTH(u);
	//IMPRIMER_POINT("u", u);
	Vector_3 v = circle.supporting_plane().base2();
	v = v / VECTOR_LENGTH(v);
	//IMPRIMER_POINT("v", v);
	double rayon = sqrt(CGAL::to_double(circle.squared_radius()));
	double pas = 2 * CGAL_PI / precision;
	glLineWidth(5.0);
	glDisable(GL_LIGHTING); // desactive la gestion de la lumiere

	Point_3 pprec = centre + rayon * u;
	Point_3 pcour;
	glBegin(GL_LINES);
	COLOR_ROUGE;
	for (double theta = pas; theta < 2 * CGAL_PI + pas; theta = theta + pas) {
		pcour = Point_3(
			centre.x() + rayon * (u.x() * cos(theta) + v.x() * sin(theta)),
			centre.y() + rayon * (u.y() * cos(theta) + v.y() * sin(theta)),
			centre.z() + rayon * (u.z() * cos(theta) + v.z() * sin(theta))
		);
		TRACE_SEGMENT(pprec, pcour);
		pprec = pcour;
	}
	glEnd();
}

}

#endif // CGAL_SPHERICAL_CIRCLE_GL__H
