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
class Qt_widget_get_circle : public Qt_widget_tool
{
public:
  typedef Point_2<R>		Point;
  typedef Circle_2<R>		Circle;
  typedef typename R::FT	FT;

  Qt_widget_get_circle() : firstpoint(false), 
			 firsttime(true){};


private:
  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == CGAL_QT_WIDGET_GET_POINT_BUTTON && !firstpoint)
    {
      FT
	x=static_cast<FT>(widget->x_real(e->x())),
	y=static_cast<FT>(widget->y_real(e->y()));
      x1 = x;
      y1 = y;
      x2 = x;
      y2 = y;
      firstpoint = TRUE;
    } else if(e->button() == CGAL_QT_WIDGET_GET_POINT_BUTTON){
      FT
	x=static_cast<FT>(widget->x_real(e->x())),
	y=static_cast<FT>(widget->y_real(e->y()));
	widget->new_object(make_object(Circle(Point(x1,y1),
		  squared_distance(Point(x1, y1), Point(x,y)))));
	firstpoint = FALSE;
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
      FT
	x=static_cast<FT>(widget->x_real(e->x())),
	y=static_cast<FT>(widget->y_real(e->y()));
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
    widget->setCursor(crossCursor);
  };

  void deactivating()
  {
    widget->setCursor(oldcursor);
    firstpoint = false;
  };

  FT	x1, //the X of the first point
		y1; //the Y of the first point
  FT	x2, //the old second point's X
		y2; //the old second point's Y
  bool	firstpoint, //true if the user left clicked once
		firsttime;  //true if the line is not drawn
};//end class 

} // namespace CGAL

#endif // CGAL_QT_WIDGET_GET_SEGMENT_H
