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
// file          : include/CGAL/IO/Qt_widget_get_polygon.h
// package       : Qt_widget
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WIDGET_GET_POLYGON_H
#define CGAL_QT_WIDGET_GET_POLYGON_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <qcursor.h>
#include <list>

namespace CGAL {

template <class Polygon>
class Qt_widget_get_polygon : public Qt_widget_layer
{
public:
  typedef typename Polygon::Point_2 Point_2;
  typedef typename Polygon::FT FT;
  
  Qt_widget_get_polygon()
    : active(false), lastx(0), lasty(0), num_points(0), poly(),
      qpoints(0), qpoints_old(0) {};

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
      FT x=static_cast<FT>(widget->x_real(e->x()));
	    FT y=static_cast<FT>(widget->y_real(e->y()));
      poly.push_back(Point_2(x,y));
      qpoints.putPoints(num_points, 1, e->x(), e->y());
      qpoints_old.putPoints(num_points, 1, e->x(), e->y());
      ++num_points;
      active=true;
      widget->setMouseTracking(TRUE);
      return;
    };
    if(e->button() == Qt::RightButton && is_pure(e->state()))
    {
      if (active)
      {
        RasterOp old_rasterop=widget->rasterOp();
        widget->painter().setRasterOp(NotROP);
        widget->painter().drawPolyline(qpoints_old);
        widget->setRasterOp(old_rasterop);
        widget->do_paint();
        FT x=static_cast<FT>(widget->x_real(e->x()));
        FT y=static_cast<FT>(widget->y_real(e->y()));
	      poly.push_back(Point_2(x, y));
        widget->new_object(make_object(poly));
		poly.erase(poly.vertices_begin(),poly.vertices_end());
        qpoints.resize(0);
        qpoints_old.resize(0);
        num_points=0;
        active=false;
      };
    };
  };

  void mouseMoveEvent(QMouseEvent *e)
  {
    if (active)
    {

      int x=(e->pos()).x();
      int y=(e->pos()).y();
      if ((lastx != x) || (lasty != y)) {
        widget->lock();
        RasterOp old_rasterop=widget->rasterOp();
        widget->painter().setRasterOp(NotROP);
        qpoints.putPoints(num_points,1,x,y);
        widget->painter().drawPolyline(qpoints);
        // Erase old polyline
        widget->painter().drawPolyline(qpoints_old);
        widget->setRasterOp(old_rasterop);
        widget->unlock();
        lastx= x;
        lasty= y;
        qpoints_old.putPoints(num_points,1,x,y);
      }
    }
  };
  void activating()
  {	
    oldcursor = widget->cursor();
    widget->setCursor(crossCursor);
  };
  
  void deactivating()
  {
    widget->setCursor(oldcursor);
  };
  
protected:
  bool active;
  int lastx, lasty;
  int num_points;
  Polygon poly;
  QPointArray qpoints;
  QPointArray qpoints_old;
private:
  QCursor oldcursor;
};

} // namespace CGAL

#endif // CGAL_QT_WIDGET_GET_POLYGON_H
