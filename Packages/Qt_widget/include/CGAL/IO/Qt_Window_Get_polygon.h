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
// file          : include/CGAL/IO/Qt_Window_Get_polygon.h
// package       : QT_window
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WINDOW_GET_POLYGON_H
#define CGAL_QT_WINDOW_GET_POLYGON_H

#include <CGAL/IO/Qt_Window.h>
#include <CGAL/IO/Qt_Window_tool.h>
#include <list>

namespace CGAL {

template <class R>
class Qt_widget_get_polygon : public Qt_widget_tool
{
  //  Q_OBJECT
public:
  typedef Point_2<R> Point;
  typedef typename R::FT FT;

  Qt_widget_get_polygon()
    : active(false), lastx(0), lasty(0), num_points(0), points(), qpoints(0), qpoints_old(0) {};

  void mousePressEvent(QMouseEvent *e)
  {
    QPainter painter(widget);
    if(e->button() == Qt::LeftButton)
      {
	FT
	  x=static_cast<FT>(widget->x_real(e->x())),
	  y=static_cast<FT>(widget->y_real(e->y()));
	points.push_back(Point(x,y));
	qpoints.putPoints(num_points,1,e->x(),e->y());
	qpoints_old.putPoints(num_points,1,e->x(),e->y());
	++num_points;
	active=true;
	widget->setMouseTracking(TRUE);
	return;
      };
    if(e->button() == Qt::RightButton)
      {
	if (active)
	  {
	    painter.drawPolyline(qpoints_old);
	    widget->update();
	    FT
	      x=static_cast<FT>(widget->x_real(e->x())),
	      y=static_cast<FT>(widget->y_real(e->y()));
	    points.push_back(Point(x,y));
	    emit new_object(make_object(points));
	    points.clear();
	    qpoints.resize(0);
	    qpoints_old.resize(0);
	    num_points=0;
	    active=false;
	    widget->setMouseTracking(FALSE);
	  };
      };
  };

  void mouseMoveEvent(QMouseEvent *e)
  {
    if (active)
      {
	QPainter painter(widget);
	int 
	  x=(e->pos()).x(),
	  y=(e->pos()).y();
	
	if ((lastx != x) || (lasty != y)) {
	  painter.setRasterOp(NotROP);
	  qpoints.putPoints(num_points,1,x,y);
	  painter.drawPolyline(qpoints);
	  // Erase old polyline
	  painter.drawPolyline(qpoints_old);
	  painter.drawPoint(300,300);
	  lastx= x;
	  lasty= y;
	  qpoints_old.putPoints(num_points,1,x,y);
	}
      }
  };

protected:
  bool active;
  int lastx, lasty;
  int num_points;
  std::list<Point> points;
  QPointArray qpoints;
  QPointArray qpoints_old;
};

} // namespace CGAL

#endif // CGAL_QT_WINDOW_GET_POLYGON_H
