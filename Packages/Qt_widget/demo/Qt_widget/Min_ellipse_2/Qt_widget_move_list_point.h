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
// file          : demo/Qt_widget/Min_ellipse_2/Qt_widget_move_list_point.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WIDGET_MOVE_LIST_POINT_H
#define CGAL_QT_WIDGET_MOVE_LIST_POINT_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>

#include <qobject.h>
#include <qpopupmenu.h>
#include <qmessagebox.h> 
#include <qcursor.h>

#include <CGAL/squared_distance_2.h> 


class Qt_widget_movepoint_helper : public CGAL::Qt_widget_layer
{
Q_OBJECT
public:
  virtual void delete_pointi(){};
  virtual void move_pointi(){};

public slots:
  void delete_point();
  void move_point();
  void stateChanged(int);
};



template <class R>
class Qt_widget_move_list_point : public Qt_widget_movepoint_helper
{
public:
  typedef typename R::Point_2	Point;
  typedef typename R::FT        FT;
  bool			        on_first,   //if the user choosed something from the popup
			        wasrepainted;//true when the widget was repainted
  Point				old_point,  //the last point stored in the list
				current_v;  //the current point
  QPopupMenu			*popup1;
  QCursor			cursor;
  std::list<Point>*		l_of_p;
  bool                          first_time;

  //constructor
  Qt_widget_move_list_point(const QCursor c=QCursor(Qt::crossCursor)) :
      on_first(FALSE), cursor(c), first_time(true) {};
  
  void pass_the_structure(std::list<Point>* l) {
    l_of_p = l;
  }
private:
  QCursor oldcursor;

  void widget_repainted(){
      wasrepainted = TRUE;
  };

  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::LeftButton && on_first)
    {
      on_first = FALSE;
    }
    if(e->button() == Qt::RightButton)
    {
      if(l_of_p->empty())
	      QMessageBox::warning( widget, "There are no points in the list!",
        "Generate some points first or add it with the \
         input tool before using this tool!");
      else{
	FT x, y;
        widget->x_real(e->x(), x);
	widget->y_real(e->y(), y);
        Point p(x, y);
        Point closest_p;  
        //this point is the closest one to the mouse coordinates
        FT min_dist;
        typename std::list<Point>::const_iterator it = l_of_p->begin();
        min_dist = CGAL::squared_distance(p, (*it));
	closest_p = (*it);
        while(it!=l_of_p->end())
        {
  	  if (min_dist > CGAL::squared_distance(p, (*it))) {
	    min_dist = CGAL::squared_distance(p, (*it));
	    closest_p = (*it);
	  }
	  it++;
	}
        RasterOp old = widget->rasterOp();	//save the initial raster mode
        widget->setRasterOp(XorROP);
        widget->lock();
          *widget << CGAL::GREEN << CGAL::PointSize (7) 
                  << CGAL::PointStyle (CGAL::DISC);
          *widget << closest_p;
        widget->unlock();
        widget->setRasterOp(old);
        popup1->popup(widget->mapToGlobal(e->pos()));
        old_point = closest_p;
        current_v = closest_p;
        wasrepainted = FALSE;
        on_first = FALSE;
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
      *widget << CGAL::GREEN << CGAL::PointSize (5) 
              << CGAL::PointStyle (CGAL::DISC);
      if(!wasrepainted)
        *widget << old_point;
      *widget << Point(x, y);
      l_of_p->remove(old_point);
      l_of_p->push_back(Point(x, y));
      widget->redraw();	//redraw the scenes
      old_point = Point(x, y);
    }
  };
  
  void activating()
  {
    oldcursor = widget->cursor();
    widget->setCursor(cursor);
    if (first_time)
    {
      popup1 = new QPopupMenu( widget, 0);
      popup1->insertItem("Delete Point", this, SLOT(delete_point()));
      popup1->insertItem("Move Point", this,  SLOT(move_point()));
      first_time = false;
    }
  };
  
  void deactivating()
  {
    widget->setCursor(oldcursor);
  };

  void delete_pointi(){
    l_of_p->remove(current_v);
    widget->redraw();	//redraw the scenes
  };
  void move_pointi(){
    on_first = TRUE;
    widget->cursor().setPos(widget->mapToGlobal(
                            QPoint(widget->x_pixel(old_point.x()), 
                            widget->y_pixel(old_point.y()))));
  };
};

#endif // CGAL_QT_WIDGET_MOVE_LIST_POINT_H
