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
// file          : include/CGAL/IO/Qt_widget_get_line.h
// package       : Qt_widget
// author(s)     : Ursu Radu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WIDGET_GET_LINE_H
#define CGAL_QT_WIDGET_GET_LINE_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_tool.h>

#ifndef CGAL_QT_WIDGET_GET_POINT_BUTTON
#define CGAL_QT_WIDGET_GET_POINT_BUTTON Qt::LeftButton
#endif


namespace CGAL {

template <class R>
class Qt_widget_get_line : public Qt_widget_tool
{
public:
  typedef Point_2<R>	Point;
  typedef Line_2<R>	Line;
  typedef typename	R::FT FT;

  Qt_widget_get_line() : firstpoint(false), 
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
      if(x1!=x || y1!=y) {
	widget->new_object(make_object(Line(Point(x1,y1),Point(x,y))));
	firstpoint = FALSE;
      }    
    }
  }

  void leaveEvent(QEvent *e)
  {
    if(firstpoint)
    {
      RasterOp old_raster = widget->rasterOp();//save the initial raster mode
      QColor old_color = widget->color();
      widget->lock();
	widget->setRasterOp(XorROP);
	*widget << CGAL::GREEN;
	*widget << Line(Point(x1,y1),Point(x2,y2));
	widget->setRasterOp(old_raster);
	widget->setColor(old_color);
      widget->unlock();
      firsttime = true;
    }
  }
  void mouseMoveEvent(QMouseEvent *e)
  {
    if(firstpoint)
    {
      FT
	x=static_cast<FT>(widget->x_real(e->x())),
	y=static_cast<FT>(widget->y_real(e->y()));
	RasterOp old_raster = widget->rasterOp();//save the initial raster mode
	QColor old_color = widget->color();
	widget->setRasterOp(XorROP);
	widget->lock();
	*widget << CGAL::GREEN;
	if(!firsttime)
	  *widget << Line(Point(x1,y1),Point(x2,y2));
	*widget << Line(Point(x1,y1),Point(x,y));
	widget->unlock();
	widget->setRasterOp(old_raster);
	widget->setColor(old_color);

	//save the last coordinates to redraw the screen
	x2 = x;
	y2 = y;
	firsttime = false;
    }
  };
  void attaching()
  {
    oldcursor = widget->cursor();
    widget->setCursor(crossCursor);
  };

  void detaching()
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

#endif // CGAL_QT_WIDGET_GET_LINE_H
