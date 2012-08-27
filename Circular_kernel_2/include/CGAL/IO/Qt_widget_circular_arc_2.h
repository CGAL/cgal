// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Radu Ursu

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_IO_QT_WIDGET_CIRCULAR_ARC_2_H
#define CGAL_IO_QT_WIDGET_CIRCULAR_ARC_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Circular_arc_2.h>
#include <cmath>

namespace CGAL {

template < typename CK >
CGAL::Qt_widget &
operator<<(CGAL::Qt_widget & widget, const CGAL::Circular_arc_2<CK> &arc)
{
    const typename CK::Circle_2 & circ = arc.supporting_circle();
    //typename CK::Circle_2 circ = arc.supporting_circle();
    const typename CK::Point_2 & center = circ.center();
    const typename CK::Circular_arc_point_2 & source = arc.source();
    const typename CK::Circular_arc_point_2 & target = arc.target();
    double rad = std::sqrt(CGAL::to_double(circ.squared_radius()));

    int x_screen   = widget.x_pixel(CGAL::to_double(center.x()));
    int y_screen   = widget.y_pixel(CGAL::to_double(center.y()));
    int x_screen_b = widget.x_pixel(CGAL::to_double(center.x()) + rad);
    int radius     = x_screen_b - x_screen;

    double a   = std::atan2( to_double(source.y() - center.y()),
		             to_double(source.x() - center.x())); 
    double a2p = std::atan2( to_double(target.y() - center.y()),
		             to_double(target.x() - center.x()));

    if (a2p <= a)
        a2p += 2 * CGAL_PI;

    double alen2 = a2p - a;

    double diff = 180/CGAL_PI*16;

    widget.get_painter().drawArc(x_screen - radius, 
				 y_screen - radius, 
				 2 * radius, 2 * radius, 
				 (int)(a * diff), 
				 (int)(alen2 * diff));
    return widget;
}

} // namespace CGAL

#endif // CGAL_IO_QT_WIDGET_CIRCULAR_ARC_2_H
