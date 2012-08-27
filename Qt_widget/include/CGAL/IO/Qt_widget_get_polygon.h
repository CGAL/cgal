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
// Author(s)     : Laurent Rineau

#ifndef CGAL_QT_WIDGET_GET_POLYGON_H
#define CGAL_QT_WIDGET_GET_POLYGON_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <qcursor.h>
#include <list>

namespace CGAL {

template <class Polygon>
class Qt_widget_get_polygon : public Qt_widget_layer
{
public:
  typedef typename Polygon::Point_2     Point_2;
  typedef typename Polygon::Segment_2   Segment_2;
  typedef typename Polygon::Edge_const_iterator
                                        ECI;
  typedef typename Polygon::Vertex_iterator
                                        VI;
  typedef typename Polygon::FT          FT;
  
  Qt_widget_get_polygon(const QCursor c=QCursor(Qt::crossCursor),
			QObject* parent = 0, const char* name = 0)
    : Qt_widget_layer(parent, name), active(false),
      first_time(true), cursor(c) {}

  void draw()
  {
    if(poly.size() > 1)
    {
      widget->lock();
      RasterOp old_rasterop=widget->rasterOp();
      widget->get_painter().setRasterOp(XorROP);
      *widget << CGAL::GREEN;
      ECI before_end = poly.edges_end();
      --before_end; // --poly.edges_end() doesn't work on g++-2.95
		    // with std::vector as the container for the polygon
      for(ECI it = poly.edges_begin(); it != before_end; it++)
        *widget << *it;
      widget->setRasterOp(old_rasterop);
      widget->unlock();
    }
    return;
  };
protected:

  bool is_pure(Qt::ButtonState s){
    if((s & Qt::ControlButton) ||
       (s & Qt::ShiftButton) ||
       (s & Qt::AltButton))
      return 0;
    else
      return 1;
  }

  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::LeftButton && is_pure(e->state()))
    {
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);

      if(!active)
      {
        active=true;
        last_of_poly = Point_2(x, y);
        poly.push_back(Point_2(x, y));	
      } else{
        if (last_of_poly == Point_2(x,y)) return;
        rubber_old = Point_2(x, y);
        if(is_simple()){
          poly.push_back(Point_2(x,y));	
          //show the last rubber as edge of the polygon
          widget->lock();
            RasterOp old_rasterop=widget->rasterOp();
            widget->get_painter().setRasterOp(XorROP);
            *widget << CGAL::WHITE;
            *widget << Segment_2(rubber, last_of_poly);
            *widget << CGAL::GREEN;
            *widget << Segment_2(rubber, last_of_poly);
            widget->setRasterOp(old_rasterop);
          widget->unlock();
          last_of_poly = Point_2(x, y);
        }
      }
      return;
    };
    if(e->button() == Qt::RightButton && is_pure(e->state()))
    {
      if (active) {
        widget->new_object(make_object(poly));
        active = false;
        first_time = true;
        poly.erase(poly.vertices_begin(), poly.vertices_end());
      }
    };
  };//end mousePressEvent

  void keyPressEvent(QKeyEvent *e)
  {
    switch ( e->key() ) {
      case Key_Escape:			// key_escape
          if(poly.size() > 1){
            widget->lock();
              RasterOp old_rasterop=widget->rasterOp();
              widget->get_painter().setRasterOp(XorROP);
              *widget << CGAL::GREEN;

	      // g++-2.95 doesn't like --poly.vertices_end() if the
	      // container of the polygon is std::vector
	      VI last_of_poly_it = poly.vertices_end();
	      --last_of_poly_it;
	      VI before_last_of_poly_it = last_of_poly_it;
	      --before_last_of_poly_it;

              *widget << Segment_2(*before_last_of_poly_it, last_of_poly);
              *widget << CGAL::WHITE;
              *widget << Segment_2(rubber, last_of_poly);
              *widget << Segment_2(rubber, *before_last_of_poly_it);
              widget->setRasterOp(old_rasterop);
            widget->unlock();
            last_of_poly = *before_last_of_poly_it; 
            poly.erase(last_of_poly_it);
          }
          break;
    }//endswitch
  }

  void mouseMoveEvent(QMouseEvent *e)
  {
    if(active) {
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      rubber = Point_2(x, y);
      if(e->state() == Qt::ShiftButton){
        FT dx, dy;
        dx = last_of_poly.x() > x ? last_of_poly.x() - x : x - last_of_poly.x();
        dy = last_of_poly.y() > y ? last_of_poly.y() - y : y - last_of_poly.y();
        widget->lock();
        RasterOp old_rasterop=widget->rasterOp();
        widget->get_painter().setRasterOp(XorROP);
        *widget << CGAL::WHITE;
        if(!first_time)
          *widget << Segment_2(rubber_old, last_of_poly);
        if(dx < dy)
          rubber = Point_2(last_of_poly.x(), y);
        else
          rubber = Point_2(x, last_of_poly.y());
        *widget << Segment_2(rubber, last_of_poly);
        widget->unlock();
        first_time = false;
        rubber_old = rubber;
        widget->cursor().setPos(widget->mapToGlobal(
          QPoint(widget->x_pixel(CGAL::to_double(rubber.x())), 
          widget->y_pixel(CGAL::to_double(rubber.y())))));
        widget->setRasterOp(old_rasterop);    
      } else { 
        widget->lock();
        RasterOp old_rasterop=widget->rasterOp();
        widget->get_painter().setRasterOp(XorROP);
        *widget << CGAL::WHITE;      	
        if(!first_time)
          *widget << Segment_2(rubber_old, last_of_poly);
        *widget << Segment_2(rubber, last_of_poly);
        first_time = false;
        rubber_old = rubber;
        widget->setRasterOp(old_rasterop);
        widget->unlock();
      }
    }
  };
  void activating()
  {
    oldcursor = widget->cursor();
    widget->setCursor(cursor);
    oldpolicy = widget->focusPolicy();
    widget->setFocusPolicy(QWidget::StrongFocus);
  };
  
  void deactivating()
  {
    poly.erase(poly.vertices_begin(), poly.vertices_end());
    active = false;
    first_time = true;
    widget->setCursor(oldcursor);
    widget->setFocusPolicy(oldpolicy);
    widget->redraw();
  };
  void leaveEvent(QEvent *)
  {
    if (active)
    {
      widget->lock();
        RasterOp old_rasterop=widget->rasterOp();
        widget->get_painter().setRasterOp(XorROP);
        *widget << CGAL::WHITE;
        *widget << Segment_2(rubber_old, last_of_poly);
        widget->setRasterOp(old_rasterop);
      widget->unlock();
      first_time = true;
    }
  }
private:
  virtual bool is_simple()
  {
    return true;
  }
  
protected:
  bool	active,     //true if the first point was inserted
        first_time; //true if it is the first time when 
                    //draw the rubber band
  Point_2 rubber,   //the new point of the rubber band
          last_of_poly, //the last point of the polygon
          rubber_old; //the old point of the rubber band
  Polygon poly;       //the polygon
  QWidget::FocusPolicy  oldpolicy;
  QCursor oldcursor;
  QCursor cursor;

};

} // namespace CGAL

#endif // CGAL_QT_WIDGET_GET_POLYGON_H
