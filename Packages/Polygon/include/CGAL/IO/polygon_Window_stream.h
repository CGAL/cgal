// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : 
// release_date  :
//
// file          : include/CGAL/IO/polygon_Window_stream.h
// source        :
// author(s)     : Wieger Wesselink, Geert-Jan Giezeman, Matthias Baesken
//
// coordinator   : Utrecht University
// ============================================================================

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
