// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Radu Ursu
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file : include/CGAL/IO/Qt_widget_circular_arc_2.h

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

    double rad = std::sqrt(CGAL::to_double(circ.squared_radius()));

    int x_screen   = widget.x_pixel(CGAL::to_double(circ.center().x()));
    int y_screen   = widget.y_pixel(CGAL::to_double(circ.center().y()));
    int x_screen_b = widget.x_pixel(CGAL::to_double(circ.center().x()) + rad);
    int radius     = x_screen_b - x_screen;

    double a   = std::atan2( to_double(arc.source().y() - circ.center().y()),
		             to_double(arc.source().x() - circ.center().x())); 
    double a2p = std::atan2( to_double(arc.target().y() - circ.center().y()),
		             to_double(arc.target().x() - circ.center().x()));

    if (a2p < a)
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
