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
// file          : snapping_layer.h
// package       : Planar_map
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : 
//
// ============================================================================

#ifndef CGAL_SNAPPING_LAYER_H
#define CGAL_SNAPPING_LAYER_H

#include "cgal_types.h"
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>

#include <qobject.h>
#include <qpopupmenu.h>
#include <qmessagebox.h> 
#include <qcursor.h>

#include <CGAL/squared_distance_2.h> 

template <class R>
class Snapping_layer : public CGAL::Qt_widget_layer
{
public:
  typedef typename R::Point_2   Point;
  typedef typename R::FT        FT;
  bool                          on_first,   //true if the user selected the first vertex
                                wasrepainted;//true when the widget was repainted
  Point                         old_point,  //the last end point found
                                current_v;  //the current point
  QCursor                       cursor;
  std::list<Curve>              *curveslist;
  Curve                         selected_curve, old_curve;
  bool                          source, target;
  Planar_map                    *pm;

  //constructor
  Snapping_layer(const QCursor c=QCursor(Qt::crossCursor)) :
      on_first(false), cursor(c),
      source(false), target(false) {};
  
  void pass_the_structure(std::list<Curve>* l, Planar_map * p) {
    curveslist = l;
    pm = p;
  }
private:
  void draw(){
      wasrepainted = true;
  };

private:
  QCursor oldcursor;

  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::LeftButton && on_first)
    {
      curveslist->remove(selected_curve);
      curveslist->push_back(old_curve);
      Locate_type i;
      Halfedge_handle hh;
      if(source)
        hh = pm->locate(selected_curve.source(), i);
      else
        hh = pm->locate(selected_curve.target(), i);
      pm->remove_edge(hh);
      pm->insert(old_curve);
      on_first = false; source = false; target = false;
      selected_curve = Curve();
      widget->redraw();
    }
    if(e->button() == Qt::RightButton && !on_first)
    {
      if(curveslist->empty())
	      QMessageBox::warning( widget, "There are no segments in the list!",
        "Input some segments before using this layer!");
      else{
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
            selected_curve = (*it); old_curve = (*it);
            source = true; target = false;
          }
          if (min_dist > CGAL::squared_distance(p, (*it).target())) {
            min_dist = CGAL::squared_distance(p, (*it).target());
            closest_p = (*it).target();
            selected_curve = (*it); old_curve = (*it);
            target = true; source = false;
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
        current_v = closest_p;
        wasrepainted = false;
        on_first = true;
      }
    }
  };

  void mouseMoveEvent(QMouseEvent *e)
  {
    if(on_first)
    {
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      RasterOp old = widget->rasterOp();	//save the initial raster mode
      widget->setRasterOp(XorROP);
      widget->lock();
      *widget << CGAL::GREEN << CGAL::PointSize (5) 
              << CGAL::PointStyle (CGAL::DISC);
      if(!wasrepainted){
        *widget << CGAL::RED << old_curve;
        *widget << CGAL::GREEN << old_point;
      }
      
      Point p(x,y);
      Point closest_p;
      FT min_dist;
      typename std::list<Curve>::const_iterator it = curveslist->begin();
      min_dist = CGAL::squared_distance(p, (*it).source());
      closest_p = (*it).source();
      
      while(it!=curveslist->end())
      {
        if( (*it) == selected_curve){
            it++;
            continue;
          }
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
      
      if(target)
        old_curve = Curve(old_curve.source(), closest_p);
      else
        old_curve = Curve(closest_p, old_curve.target());
      *widget << CGAL::RED << old_curve;
      *widget << CGAL::GREEN << closest_p;
      widget->unlock();
      widget->setRasterOp(old);
      old_point = closest_p;
    }
  }; 

  void activating()
  {
    oldcursor = widget->cursor();
    widget->setCursor(cursor);
    on_first = false; source = false;
    target = false;
  };
  
  void deactivating()
  {
    widget->setCursor(oldcursor);
  };

};

#endif // CGAL_SNAPPING_LAYER_H
