// ============================================================================
//
// Copyright (c) 1997-2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : segment_input_layer_with_snapping.h
// package       : Planar_map
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : 
//
// ============================================================================


#include "cgal_types.h"
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_get_segment.h>



template <class R>
class Segment_input_layer : public CGAL::Qt_widget_get_segment<R>{
private:
  std::list<Curve>              *curveslist;
  Planar_map                    *pm;
  //true if the user selected the first vertex
  Point                         old_point;
  QCursor                       cursor;
public:
  typedef typename R::Segment_2	Segment;
  typedef typename R::Point_2   Point;
  typedef typename R::FT        FT;
public:
 //constructor
  Segment_input_layer(const QCursor c=QCursor(Qt::crossCursor)) :
        CGAL::Qt_widget_get_segment<R>(), cursor(c){}
  void pass_the_structure(std::list<Curve>* l, Planar_map * p) {
    curveslist = l;
    pm = p;
  }
private:
  void mousePressEvent(QMouseEvent *e)
  {
    CGAL::Qt_widget_get_segment<Kernel>::mousePressEvent(e);
    if(e->button() == Qt::RightButton && is_pure(e->state()))
    {
      if(curveslist->empty()){
	      QMessageBox::warning( widget, "There are no segments in the list!",
        "Input some segments using the left mouse button!");
        return;
      }
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      Point p(x, y);
      Point closest_p;  
      //this point is the closest one to the mouse coordinates
      FT min_dist;
      typename std::list<Curve>::const_iterator it = curveslist->begin();
      min_dist = CGAL::squared_distance(p, (*it).source());
      closest_p = (*it).source();
      
      while(it!=curveslist->end())
      {
        if (min_dist > CGAL::squared_distance(p, (*it).source())) {
          min_dist = CGAL::squared_distance(p, (*it).source());
          closest_p = (*it).source();
        }
        if (min_dist > CGAL::squared_distance(p, (*it).target())) {
          min_dist = CGAL::squared_distance(p, (*it).target());
          closest_p = (*it).target();
        }
        it++;
      }
      
      RasterOp old = widget->rasterOp();	//save the initial raster mode
      widget->setRasterOp(XorROP);
      widget->lock();
        *widget << CGAL::GREEN << CGAL::PointSize (5) 
                << CGAL::PointStyle (CGAL::DISC);
        *widget << closest_p;
      widget->unlock();
      widget->setRasterOp(old);        
      old_point = closest_p;
      if(!firstpoint){
        x1 = closest_p.x();
        y1 = closest_p.y();
        x2 = closest_p.x();
        y2 = closest_p.y();
        firstpoint = true;
      } else {
        if(x1 != closest_p.x() || y1 != closest_p.y()) {
          widget->new_object(
            CGAL::make_object(Segment(Point(x1, y1), 
              Point(closest_p.x(), closest_p.y()))));
          firstpoint = false;
        }
      }
    }
  }
  void activating()
  {
    oldcursor = widget->cursor();
    widget->setCursor(cursor);
    firstpoint = false;
  };
  
  void deactivating()
  {
    widget->setCursor(oldcursor);
  };
};
