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
// file          : include/CGAL/IO/Qt_widget_MovePoint.h
// package       : QT_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


#ifndef CGAL_QT_WIDGET_MOVEPOINT_H
#define CGAL_QT_WIDGET_MOVEPOINT_H


#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>


#include <qobject.h>
#include <qpopupmenu.h>
#include <qcursor.h>

#include <list>

class Qt_widget_movepoint_helper : public CGAL::Qt_widget_layer
{
Q_OBJECT
public:
virtual void delete_pointi(){};
virtual void move_pointi(){};

public slots:
  void stateChanged(int i);
  void delete_point();
  void move_point();
};

template <class T, class A>
class Qt_widget_movepoint : public Qt_widget_movepoint_helper
{
public:

  typedef typename T::Point                       Point;
  typedef typename T::Segment			  Segment;
  typedef typename T::Face_handle                 Face_handle;
  typedef typename T::Vertex_handle               Vertex_handle;
  typedef typename T::Vertex_iterator             Vertex_iterator;
  typedef typename T::Geom_traits::FT             FT;
  typedef std::list<Point>                        CGALPointlist;
protected:
  FT            first_x, first_y;
  FT            x2, y2;
  bool          wasrepainted;
  bool          on_first;

  Vertex_handle current_v;	//the vertex that will be process
  Point	  old_point;
  T             *dt;
  A             *as;
  QPopupMenu    *popup1;
  CGALPointlist L;

public:
  Qt_widget_movepoint() : wasrepainted(true), on_first(FALSE)
  {};
  void set_variables (T *t, A* a) {dt = t; as = a;}
private:
  QCursor oldcursor;

  void draw(){
    wasrepainted = true;
  };
  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::LeftButton && on_first)
    {
      on_first = FALSE;
    }
    if(e->button() == Qt::RightButton)
    {
      if (dt->dimension()<2) return;
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      Point p(x, y);
      Vertex_handle v = dt->nearest_vertex(p);
      RasterOp old = widget->rasterOp();	//save the initial raster mode
      widget->setRasterOp(XorROP);
      widget->lock();
        *widget << CGAL::GREEN << CGAL::PointSize (7) 
              << CGAL::PointStyle (CGAL::DISC);
	      if(!wasrepainted)
          *widget << old_point;
          *widget << v->point();
      widget->unlock();
      widget->setRasterOp(old);
      popup1->popup(widget->mapToGlobal(e->pos()));
      old_point = v->point();
      current_v = v;
      wasrepainted = FALSE;
      on_first = FALSE;
    }	
  };
  void mouseMoveEvent(QMouseEvent *e)
  {
    if(on_first)
    {
      FT x, y;
      widget->x_real(e->x(), x),
      widget->y_real(e->y(), y);
  		
      *widget << CGAL::GREEN << CGAL::PointSize (5)
              << CGAL::PointStyle (CGAL::DISC);
      if(!wasrepainted)
        *widget << old_point;
      *widget << Point(x, y);
      dt->remove(current_v);
      current_v = dt->insert(Point(x, y));
      FT alpha_index = as->get_alpha();
      as->clear();
      L.clear();
      Vertex_iterator it = dt->vertices_begin(), 
	              beyond = dt->vertices_end();
      while(it != beyond) {      
        L.push_back((*it).point());
        ++it;
      }
      as->make_alpha_shape(L.begin(), L.end());
      as->set_alpha(alpha_index);
      widget->redraw();	//redraw the layers
      old_point = Point(x, y);
    }
  }
  void activating()
  {
    oldcursor = widget->cursor();
    widget->setCursor(crossCursor);
    wasrepainted = TRUE;
    popup1 = new QPopupMenu( widget, 0);
    popup1->insertItem("Delete Point", this, SLOT(delete_point()));
    popup1->insertItem("Move Point", this,  SLOT(move_point()));
  }

  void deactivating()
  {
    widget->setCursor(oldcursor);
  }
  void delete_pointi(){
    dt->remove(current_v);
    FT alpha_index = as->get_alpha();
    as->clear();
    L.clear();
    Vertex_iterator it = dt->vertices_begin(), 
	            beyond = dt->vertices_end();
    while(it != beyond) {      
      L.push_back((*it).point());
      ++it;
    }
    as->make_alpha_shape(L.begin(), L.end());
    as->set_alpha(alpha_index);
    widget->redraw();	//redraw the scenes
  }
  void move_pointi(){
    on_first = TRUE;
    widget->cursor().setPos(widget->mapToGlobal(
                            QPoint(widget->x_pixel(old_point.x()), 
                            widget->y_pixel(old_point.y()))));
  }
};//end class 


#endif // CGAL_QT_WIDGET_GET_SEGMENT_H

