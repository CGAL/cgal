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
// file          : include/CGAL/IO/Qt_widget_get_circle.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WIDGET_GET_CIRCLE_H
#define CGAL_QT_WIDGET_GET_CIRCLE_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/squared_distance_2.h> 

#ifndef CGAL_QT_WIDGET_GET_POINT_BUTTON
#define CGAL_QT_WIDGET_GET_POINT_BUTTON Qt::LeftButton
#endif


namespace CGAL {

template <class R>
class Qt_widget_get_circle : public Qt_widget_layer
{
public:
  typedef typename R::Point_2		Point;
  typedef typename R::Circle_2		Circle;
  typedef typename R::FT	FT;

  Qt_widget_get_circle(const QCursor c=QCursor(Qt::crossCursor),
		       QObject* parent = 0, const char* name = 0)
     : Qt_widget_layer(parent, name), cursor(c), firstpoint(false),
       firsttime(true){};

  void draw(){
    firsttime = true;
  }

private:
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
    if(e->button() == CGAL_QT_WIDGET_GET_POINT_BUTTON 
       && !firstpoint
       && is_pure(e->state()))
    {
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      x1 = x;
      y1 = y;
      x2 = x;
      y2 = y;
      firstpoint = true;
    } else if(e->button() == CGAL_QT_WIDGET_GET_POINT_BUTTON){
      FT x, y; 
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      widget->new_object(make_object(Circle(Point(x1,y1),
	      squared_distance(Point(x1, y1), Point(x,y)))));
      firstpoint = false;
      firsttime = true;
    }
  };

  void leaveEvent(QEvent *)
  {
    if(firstpoint)
    {
      QColor old_color = widget->color();
      RasterOp old_raster = widget->rasterOp();//save the initial raster mode
      
      widget->lock();
        widget->setRasterOp(XorROP);
	*widget << CGAL::GREEN;
        *widget << Circle(Point(x1,y1),
        squared_distance(Point(x1, y1), Point(x2,y2)));
      widget->unlock();
      widget->setRasterOp(old_raster);
      widget->setColor(old_color);
      firsttime = true;
    }
  }

  void mouseMoveEvent(QMouseEvent *e)
  {
    if(firstpoint==TRUE)
    {		
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      QColor old_color = widget->color();
      RasterOp old_raster = widget->rasterOp();//save the initial raster mode		
      widget->setRasterOp(XorROP);
      widget->lock();
      *widget << CGAL::GREEN;
      if(!firsttime)
	*widget << Circle(Point(x1,y1),
		  squared_distance(Point(x1, y1), Point(x2,y2)));
      *widget << Circle(Point(x1,y1),
		  squared_distance(Point(x1, y1), Point(x,y)));
      widget->unlock();
      widget->setRasterOp(old_raster);
      widget->setColor(old_color);

      //save the last coordinates to redraw the screen
      x2 = x;
      y2 = y;	
      firsttime = false;
    }
  };
  void activating()
  {
    oldcursor = widget->cursor();
    widget->setCursor(cursor);
  };
  
  void deactivating()
  {
    widget->setCursor(oldcursor);
    firstpoint = false;
  };

  QCursor cursor;
  QCursor oldcursor;
  

  FT	x1, //the X of the first point
		  y1; //the Y of the first point
  FT	x2, //the old second point's X
		  y2; //the old second point's Y
  bool	firstpoint, //true if the user left clicked once
	firsttime;  //true if the line is not drawn
};//end class 

} // namespace CGAL

#endif // CGAL_QT_WIDGET_GET_SEGMENT_H
