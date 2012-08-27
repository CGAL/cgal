// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Laurent Rineau and Radu Ursu

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
  typedef Qt_widget_get_polygon<Polygon>  Get_polygon;
  typedef typename Polygon::Point_2       Point_2;
  typedef typename Polygon::Segment_2     Segment_2;
  typedef typename Polygon::Edge_const_iterator  ECI;

  Qt_widget_get_simple_polygon(const QCursor
                               c=QCursor(Qt::crossCursor),QObject*
                               parent = 0, const char* name = 0)
    : Qt_widget_get_polygon<Polygon>(c, parent, name){}
  
protected:

  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::RightButton && this->is_pure(e->state()))
    {
      if (this->active) {
        if(!this->poly.is_simple()) return;
        if(this->poly.is_clockwise_oriented())
          this->poly.reverse_orientation ();
        CGAL_assertion( ! this->poly.is_clockwise_oriented());
      }
    }
    Get_polygon::mousePressEvent(e);
  }; 

private:
  bool is_simple()
  {
    Segment_2 rubber_segment(this->rubber, this->last_of_poly);
    if(this->poly.size() > 1)
    {
      ECI before_last_it = this->poly.edges_end();
      --before_last_it;
      --before_last_it;
      ECI it;
      for(it = this->poly.edges_begin(); it != before_last_it; it++)
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
