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
// file          : demo/Qt_widget/Triangulation_2/triangulation_2_layers.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_2_LAYERS_H
#define CGAL_TRIANGULATION_2_LAYERS_H

#include <CGAL/Cartesian.h>
#include <CGAL/Circle_2.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_Triangulation_2.h>
#include <qobject.h>


template <class T>
class Qt_layer_show_triangulation : public CGAL::Qt_widget_layer
{
public:
	
  Qt_layer_show_triangulation(T &t) : tr(t){};


  void draw()
  {
    *widget << CGAL::BLUE; 
    *widget << tr;
  };
	
private:
  T &tr;
};//end class 

template <class T>
class Qt_layer_show_voronoi : public CGAL::Qt_widget_layer
{
public:
  Qt_layer_show_voronoi(T &t1) : tr(t1){};

  void draw()
  {
    *widget << CGAL::RED ;
    tr.draw_dual(*widget);
  };
	
private:
  T	&tr;
};//end class 

template <class T>
class Qt_layer_show_points : public CGAL::Qt_widget_layer {
public:
  typedef typename T::Point           Point;
  typedef typename T::Segment         Segment;
  typedef typename T::Vertex          Vertex;
  typedef typename T::Vertex_iterator	Vertex_iterator;

  Qt_layer_show_points(T &t) : tr(t){};

  void draw()
  {  
    Vertex_iterator it = tr.vertices_begin(), 
		beyond = tr.vertices_end();
    *widget << CGAL::GREEN << CGAL::PointSize (3) 
		<< CGAL::PointStyle (CGAL::DISC);    
    while(it != beyond) {      
      *widget << (*it).point();
      ++it;
    }
  };
private:
  T	&tr;
  
};//end class 

template <class T>
class Qt_layer_nearest_vertex : public CGAL::Qt_widget_layer
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
  void leaveEvent(QEvent *)
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


template <class T>
class Qt_layer_circum_circle : public CGAL::Qt_widget_layer
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
    FT x, y;
    widget->x_real(e->x(), x);
    widget->y_real(e->y(), y);
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
  void leaveEvent(QEvent *)
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


#endif
