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
// file          : include/CGAL/IO/Qt_layer_show_nearest_vertex.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_LAYER_SHOW_NEAREST_VERTEX_H
#define CGAL_QT_LAYER_SHOW_NEAREST_VERTEX_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <qobject.h>


namespace CGAL {

template <class T>
class Qt_layer_nearest_vertex : public Qt_widget_layer
{
public:
  typedef typename T::Point			Point;
  typedef typename T::Segment			Segment;
  typedef typename T::Face_handle		Face_handle;
  typedef typename T::Vertex_handle		Vertex_handle;
  typedef typename T::Geom_traits::FT		FT;

  Qt_layer_nearest_vertex(T &t) : tr(t), first_time(TRUE) {};
	
  void draw(){first_time = TRUE;};
  void mouseMoveEvent(QMouseEvent *e)
  {
    if (tr.dimension()<1) return;
    FT x, y;
    widget->x_real(e->x(), x),
    widget->y_real(e->y(), y);
    Point p(x, y);
    RasterOp old = widget->rasterOp();	//save the initial raster mode
    widget->setRasterOp(XorROP);
    widget->lock();
    Vertex_handle v = tr.nearest_vertex(p);
    *widget << CGAL::GREEN << CGAL::PointSize (10)
		<< CGAL::PointStyle (CGAL::CIRCLE);
    if(!first_time)
      *widget << oldPoint;	
    *widget << v->point();
    widget->unlock();
    widget->setRasterOp(old);
    oldPoint = v->point();
    first_time = FALSE;
  };
  void leaveEvent(QEvent *e)
  {
    widget->lock();
    RasterOp old = widget->rasterOp();	//save the initial raster mode
    widget->setRasterOp(XorROP);
    *widget << CGAL::GREEN << CGAL::PointSize (10) 
		<< CGAL::PointStyle (CGAL::CIRCLE);
    *widget << oldPoint;	
    widget->unlock();
    widget->setRasterOp(old);
    first_time = TRUE;
    //remove_leftovers(widget);
  }

private:
	T     &tr;
	Point oldPoint, newPoint;
	bool  first_time;
	
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_NEAREST_VERTEX_H
