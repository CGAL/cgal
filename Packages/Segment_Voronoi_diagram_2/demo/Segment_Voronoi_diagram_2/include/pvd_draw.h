// Copyright (c) 2004,2005  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>


#ifndef SVD_DRAW_H
#define SVD_DRAW_H



template<class T, class Widget>
void draw_diagram(Widget& widget, const T& svd)
{
  widget << CGAL::BLUE;
#if !defined (__POWERPC__)
  widget << CGAL::PointSize(3);
  widget << CGAL::LineWidth(3);
#endif

  typename T::Finite_edges_iterator eit = svd.finite_edges_begin();
  for (; eit != svd.finite_edges_end(); ++eit) {
    if ( eit->first->vertex( svd.cw(eit->second) )->info() !=
	 eit->first->vertex( svd.ccw(eit->second) )->info() ) {
      svd.draw_dual_edge(*eit, widget);	
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
      svd.draw_dual_edge(*eit, widget);	
    }
#endif
  }
}




#endif // SVD_DRAW_H
