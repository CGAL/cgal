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
#include <CGAL/IO/Qt_widget_tool.h>
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

  Qt_widget_get_circle() {};

private:
  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == CGAL_QT_WIDGET_GET_POINT_BUTTON)
    {
      FT
	x=static_cast<FT>(widget->x_real(e->x())),
	y=static_cast<FT>(widget->y_real(e->y()));
      x1 = x;
      y1 = y;
      x2 = x;
      y2 = y;
      firstpoint = TRUE;
    }
  };

  void mouseMoveEvent(QMouseEvent *e)
  {
    if(firstpoint==TRUE)
    {		
      FT
	x=static_cast<FT>(widget->x_real(e->x())),
	y=static_cast<FT>(widget->y_real(e->y()));
	  QColor oldcolor = widget->color();
      RasterOp old = widget->rasterOp();	//save the initial raster mode		
      widget->setRasterOp(XorROP);
      widget->lock();
      *widget << CGAL::GREEN;
      *widget << Circle(Point(x1,y1),
						squared_distance(Point(x1, y1), Point(x2,y2)));
      *widget << Circle(Point(x1,y1),
						squared_distance(Point(x1, y1), Point(x,y)));
      widget->unlock();
      widget->setRasterOp(old);
	  widget->setColor(oldcolor);

      //save the last coordinates to redraw the screen
      x2 = x;
      y2 = y;	
    }
  };
  void mouseReleaseEvent(QMouseEvent *e)
  {
    if(e->button() == CGAL_QT_WIDGET_GET_POINT_BUTTON)		
    {
      FT
	x=static_cast<FT>(widget->x_real(e->x())),
	y=static_cast<FT>(widget->y_real(e->y()));
    //  if(x1!=x || y1!=y)
	widget->new_object(make_object(Circle(Point(x1,y1),
				squared_distance(Point(x1, y1), Point(x,y)))));
	firstpoint = FALSE;
    }
  }

  void attaching()
  {
    oldcursor = widget->cursor();
    widget->setCursor(crossCursor);
  };

  void detaching()
  {
    widget->setCursor(oldcursor);
  
  };

  FT	x1, y1;
  FT	x2, y2;
  bool	firstpoint;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_WIDGET_GET_SEGMENT_H
