// Copyright (c) 2004,2005  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>


#ifndef PDG_DRAW_H
#define PDG_DRAW_H



template<class T, class Widget>
void draw_diagram(Widget& widget, const T& sdg)
{
  widget << CGAL::BLUE;
#if !defined (__POWERPC__)
  widget << CGAL::PointSize(3);
  widget << CGAL::LineWidth(3);
#endif

  typename T::Finite_edges_iterator eit = sdg.finite_edges_begin();
  for (; eit != sdg.finite_edges_end(); ++eit) {
    if ( eit->first->vertex( sdg.cw(eit->second) )->info() !=
	 eit->first->vertex( sdg.ccw(eit->second) )->info() ) {
      sdg.draw_dual_edge(*eit, widget);
    }
#if 0
    Site_2 p = eit->first->vertex(  cw(eit->second) )->site();
    Site_2 q = eit->first->vertex( ccw(eit->second) )->site();

    bool is_endpoint_of_seg =
      ( p.is_segment() && q.is_point() &&
	is_endpoint_of_segment(q, p) ) ||
      ( p.is_point() && q.is_segment() &&
	is_endpoint_of_segment(p, q) );

    if ( !is_endpoint_of_seg ) {
      sdg.draw_dual_edge(*eit, widget);
    }
#endif
  }
}




#endif // PDG_DRAW_H
