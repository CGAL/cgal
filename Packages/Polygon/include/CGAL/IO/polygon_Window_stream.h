// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Wieger Wesselink, Geert-Jan Giezeman, Matthias Baesken

// Window_stream I/O operators
// ===========================

// Polygon_2
// ---------
#ifdef CGAL_POLYGON_2_H
#ifndef CGAL_WINDOW_STREAM_POLYGON_2_H
#define CGAL_WINDOW_STREAM_POLYGON_2_H


#if defined(CGAL_USE_CGAL_WINDOW)
#include <list>
#endif


CGAL_BEGIN_NAMESPACE

template <class Traits, class Container>
Window_stream&
operator<<(Window_stream& ws,
           const Polygon_2<Traits,Container> &polygon)
{
  typedef typename Polygon_2<Traits,Container>::Edge_const_circulator EI;
  typedef typename Traits::Point_2  Point_2;
  
  EI e = polygon.edges_circulator();
  
#if defined(CGAL_USE_CGAL_WINDOW)
  CGAL::color cl = ws.get_fill_color();
  
  if (cl != CGAL::invisible) { 
    std::list<CGAL::window_point> LP;
      
    if (e != NULL) {
     EI end = e;
     do {
      Point_2 p = (*e).source();
      double x = CGAL::to_double(p.x());
      double y = CGAL::to_double(p.y());      
     
      LP.push_back(CGAL::window_point(x,y));
      ++e;
      } while (e != end);
    }
    
    ws.draw_filled_polygon(LP,cl);    
  }
#else  
  leda_color cl = ws.get_fill_color();
  
  if (cl != leda_invisible) { // draw filled polygon ...
    leda_list<leda_point> LP;
      
    if (e != NULL) {
     EI end = e;
     do {
      Point_2 p = (*e).source();
      double x = CGAL::to_double(p.x());
      double y = CGAL::to_double(p.y());      
     
      LP.append(leda_point(x,y));
      ++e;
      } while (e != end);
    }
    
    ws.draw_filled_polygon(LP,cl);    
  }
#endif  
  else {
  if (e != NULL) {
    EI end = e;
    do {
      ws << *e;
      ws << (*e).source();
      ++e;
      } while (e != end);
    }
  }
  
  return ws;
}

CGAL_END_NAMESPACE

#endif // CGAL_WINDOW_STREAM_POLYGON_2_H
#endif // CGAL_POLYGON_2_H
