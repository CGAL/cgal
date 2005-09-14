// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file : include/CGAL/IO/Qt_widget_circular_arc_endpoint_2.h

#ifndef CGAL_IO_QT_WIDGET_CIRCULAR_ARC_ENDPOINT_2_H
#define CGAL_IO_QT_WIDGET_CIRCULAR_ARC_ENDPOINT_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Circular_arc_endpoint_2.h>

namespace CGAL {

template < typename CK >
CGAL::Qt_widget &
operator<<(CGAL::Qt_widget & widget, const CGAL::Circular_arc_point_2<CK> &p)
{
  typedef typename CK::Point_2   Point_2;
  return widget << Point_2(to_double(p.x()), to_double(p.y()));
}

} // namespace CGAL

#endif // CGAL_IO_QT_WIDGET_CIRCULAR_ARC_ENDPOINT_2_H
