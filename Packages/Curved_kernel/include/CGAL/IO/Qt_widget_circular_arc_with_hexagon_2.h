// Copyright (c) 2005  INRIA Sophia-Antipolis (France) and
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

#ifndef CGAL_IO_QT_WIDGET_HEX_CIRCULAR_ARC_2_H
#define CGAL_IO_QT_WIDGET_HEX_CIRCULAR_ARC_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Filtered_hexagon_curved_kernel/Circular_arc_with_hexagon_2.h>
#include <cmath>

namespace CGAL {

template < typename CK >
Qt_widget &
operator<<(Qt_widget & widget, 
           const Circular_arc_with_hexagon_2<CK> &a)
{    
  typedef Circular_arc_with_hexagon_2<CK>                       Circular_arc_2;
  typedef typename Circular_arc_2::Hexagon                      Hexagon;

  widget<<CGAL::ORANGE;
  widget<<a.arc();

  typename Circular_arc_2::Hexagon_const_iterator hix1;

  for (hix1=a.hexagons_begin(); hix1!= a.hexagons_end(); ++hix1)
    {
      widget<<CGAL::BLACK;

      typename Hexagon::Edge_const_iterator it1;
      typename Hexagon::Vertex_const_iterator vt1;

      for(it1=(*hix1).edges_begin(); it1!=(*hix1).edges_end(); ++it1)
        widget<<*it1;

      widget<<CGAL::GREEN;

      for(vt1=(*hix1).vertices_begin(); vt1!=(*hix1).vertices_end(); ++vt1)
        widget<<*vt1;

    }

    return widget;
}

} // namespace CGAL

#endif // CGAL_IO_QT_WIDGET_HEX_CIRCULAR_ARC_2_H
