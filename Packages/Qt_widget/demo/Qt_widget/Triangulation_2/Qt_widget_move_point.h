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

namespace CGAL {
  class Qt_widget_movepoint_helper : public Qt_widget_layer
  {
	Q_OBJECT
  public:
	virtual void delete_pointi(){};
	virtual void move_pointi(){};

  public slots:
    void delete_point();
    void move_point();
  };

  template <class T>
  class Qt_widget_movepoint : public Qt_widget_movepoint_helper
  {
  public:

    typedef typename T::Point			      Point;
    typedef typename T::Segment			    Segment;
    typedef typename T::Face_handle		  Face_handle;
    typedef typename T::Vertex_handle		Vertex_handle;
    typedef typename T::Geom_traits::FT	FT;
  protected:
    FT						first_x, first_y;
    FT						x2, y2;
    bool					wasrepainted;
    bool					on_first;

    Vertex_handle			current_v;	//the vertex that will be process
    Point					    old_point;
    T						      *dt;
    QPopupMenu				*popup1;
  public:
    Qt_widget_movepoint() : wasrepainted(true), on_first(FALSE)
    {
      
    };
    void set_Delaunay (T *t) {dt = t;}
  private:
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
	        FT
	          x=static_cast<FT>(widget->x_real(e->x())),
	          y=static_cast<FT>(widget->y_real(e->y()));

	      Point p(x, y);
	      Vertex_handle v = dt->nearest_vertex(p);
	      RasterOp old = widget->rasterOp();	//save the initial raster mode
	      widget->setRasterOp(XorROP);
	      widget->lock();
	      *widget << CGAL::GREEN << CGAL::PointSize (7) << CGAL::PointStyle (CGAL::DISC);
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
      FT
	      x=static_cast<FT>(widget->x_real(e->x())),
	      y=static_cast<FT>(widget->y_real(e->y()));
			
      *widget << CGAL::GREEN << CGAL::PointSize (5) << CGAL::PointStyle (CGAL::DISC);
      if(!wasrepainted)
        *widget << old_point;
      *widget << Point(x, y);
      dt->remove(current_v);
      current_v = dt->insert(Point(x, y));
      widget->redraw();	//redraw the scenes
      old_point = Point(x, y);
    }
  };
    void activating()
    {
      oldcursor = widget->cursor();
      widget->setCursor(crossCursor);
      wasrepainted = TRUE;
      popup1 = new QPopupMenu( widget, 0);
      popup1->insertItem("Delete Point", this, SLOT(delete_point()));
      popup1->insertItem("Move Point", this,  SLOT(move_point()));
    };

    void deactivating()
    {
      widget->setCursor(oldcursor);
    };
    void delete_pointi(){
      dt->remove(current_v);
      widget->redraw();	//redraw the scenes
    };
    void move_pointi(){
      on_first = TRUE;
    };
  };//end class 



} // namespace CGAL



#endif // CGAL_QT_WIDGET_GET_SEGMENT_H

