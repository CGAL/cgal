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
// file          : include/CGAL/IO/Qt_layer_show_circum_circle.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_LAYER_SHOW_CIRCUM_CIRCLE_H
#define CGAL_QT_LAYER_SHOW_CIRCUM_CIRCLE_H

#include <CGAL/Cartesian.h>
#include <CGAL/Circle_2.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <qobject.h>


namespace CGAL {

template <class T>
class Qt_layer_circum_circle : public Qt_widget_layer
{
public:
  typedef typename T::Point			Point;
  typedef typename T::Segment			Segment;
  typedef typename T::Face_handle		Face_handle;
  typedef typename T::Vertex_handle		Vertex_handle;
  typedef typename T::Geom_traits::FT		FT;
  typedef typename CGAL::Cartesian<FT> Rep;
  typedef typename CGAL::Circle_2<Rep>           Circle;

  Qt_layer_circum_circle(T &t) : tr(t), first_time(TRUE) {};
	
  void draw(){first_time = TRUE;};
  void mouseMoveEvent(QMouseEvent *e)
  {
    if (tr.dimension()<2) return;
    FT x=static_cast<FT>(widget->x_real(e->x()));
    FT y=static_cast<FT>(widget->y_real(e->y()));
    Point p(x, y); //that's where the mouse is situated
    Face_handle f = tr.locate(p);
    if(!tr.is_infinite(f)){
      Point circum_center = tr.dual(f);
      Vertex_handle v1 = (*f).vertex(1);
      newCircle = Circle(circum_center, 
			 CGAL::squared_distance((*v1).point(), 
						circum_center));
      RasterOp old = widget->rasterOp();	//save the initial raster mode
      widget->setRasterOp(XorROP);
      widget->lock();
      *widget << CGAL::GREEN << CGAL::PointSize(10);
      if(!first_time)
	*widget << oldCircle;	
      *widget << newCircle;
      widget->unlock();
      widget->setRasterOp(old);
      oldCircle = newCircle;
      first_time = FALSE;
    } else if(!first_time){
      RasterOp old = widget->rasterOp();	//save the initial raster mode
      widget->setRasterOp(XorROP);
      widget->lock();
      *widget << CGAL::GREEN << CGAL::PointSize(10);
      *widget << oldCircle;
      widget->unlock();
      widget->setRasterOp(old);
      first_time = true;
    }
  };
  void leaveEvent(QEvent *e)
  {
    widget->lock();
    RasterOp old = widget->rasterOp();	//save the initial raster mode
    widget->setRasterOp(XorROP);
    *widget << CGAL::GREEN;
    if(!first_time)
      *widget << oldCircle;	
    widget->unlock();
    widget->setRasterOp(old);
    first_time = TRUE;
    //remove_leftovers(widget);
  }

private:
	T       &tr;
	Circle oldCircle, newCircle;
	bool		first_time;
	
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_CIRCUM_CIRCLE_H
