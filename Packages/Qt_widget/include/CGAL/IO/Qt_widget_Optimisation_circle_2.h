// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : include/CGAL/IO/Qt_widget_Optimisation_circle_2.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WIDGET_OPTIMISATION_CIRCLE_2_H
#define CGAL_QT_WIDGET_OPTIMISATION_CIRCLE_2_H

#include <CGAL/IO/Qt_widget.h>

namespace CGAL{
template<class Traits>
Qt_widget&
operator << (Qt_widget &ws,
	     const CGAL::Optimisation_circle_2<Traits>& oc){
  typedef typename Traits::Point_2  Point_2;
  typedef typename Traits::Circle_2 Circle_2;

  double cx( CGAL::to_double( oc.center().x()));
  double cy( CGAL::to_double( oc.center().y()));
  double sr( CGAL::to_double( oc.squared_radius()));

  if( ! CGAL_NTS is_negative(sr))
    ws << Circle_2( Point_2(cx, cy), sr);
  return ws;
}

}//end namespace

#endif
