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
// file          : include/CGAL/IO/Qt_widget_get_simple_polygon.h
// package       : Qt_widget
// author(s)     : Laurent Rineau && Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WIDGET_GET_SIMPLE_POLYGON_H
#define CGAL_QT_WIDGET_GET_SIMPLE_POLYGON_H

#include <CGAL/IO/Qt_widget_get_polygon.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>  
#include <list>

#include <qcursor.h>

namespace CGAL {
template <class Polygon>
class Qt_widget_get_simple_polygon : public Qt_widget_get_polygon<Polygon>
{
public:
  typedef Qt_widget_get_polygon<Polygon> Get_polygon;

  typedef typename Polygon::Point_2   Point_2;
  typedef typename Polygon::Segment_2 Segment_2;
  typedef typename Polygon::Edge_const_iterator  ECI;
  
  Qt_widget_get_simple_polygon() {}

protected:

  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::RightButton && is_pure(e->state()))
    {
      if (active) {
        if(!poly.is_simple()) return;
        if(poly.is_clockwise_oriented())
          poly.reverse_orientation ();
        assert( ! poly.is_clockwise_oriented());
      }
    }
    Get_polygon::mousePressEvent(e);
  }; 

private:
  bool is_simple()
  {
    Segment_2 rubber_segment(rubber, last_of_poly);
    if(poly.size() > 1)
    {
      ECI before_last_it = poly.edges_end();
      --before_last_it;
      --before_last_it;
      ECI it;
      for(it = poly.edges_begin(); it != before_last_it; it++)
      {
        if(do_intersect(*it, rubber_segment))
	  return false;
      }
      //if I'm out of this means that all the edges, 
      //didn't intersect the last one
      ++it;
      Object o = intersection(*it, rubber_segment);
      Point_2 p;
      if(assign(p, o))
	return true;
      else
	return false;
    }
    else
      return true;
  }
};

} // namespace CGAL

#endif // CGAL_QT_WIDGET_GET_SIMPLE_POLYGON_H
