// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Radu Ursu

#ifndef CGAL_QT_WIDGET_GET_SEGMENT_H
#define CGAL_QT_WIDGET_GET_SEGMENT_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>

#ifndef CGAL_QT_WIDGET_GET_POINT_BUTTON
#define CGAL_QT_WIDGET_GET_POINT_BUTTON Qt::LeftButton
#endif


namespace CGAL {

template <class R>
class Qt_widget_get_segment : public Qt_widget_layer
{
public:
  typedef typename R::Point_2		Point;
  typedef typename R::Segment_2		Segment;
  typedef typename R::FT	FT;

  Qt_widget_get_segment(const QCursor c=QCursor(Qt::crossCursor),
			QObject* parent = 0, const char* name = 0)
    : Qt_widget_layer(parent, name), cursor(c), firstpoint(false),
      firsttime(true){};

protected:
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
      firstpoint = TRUE;
    } else if(e->button() == CGAL_QT_WIDGET_GET_POINT_BUTTON 
              && is_pure(e->state())){
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      if(x1!=x || y1!=y) {
        widget->new_object(
          make_object(Segment(Point(x1,y1),Point(x,y))));
        firstpoint = FALSE;
      }    
    }
  }

  void keyPressEvent(QKeyEvent *e)
  {
    switch ( e->key() ) {
      case Key_Escape:			// key_escape
         if(firstpoint)
         {
	   firstpoint = false;
           RasterOp old_raster = widget->rasterOp();
           QColor old_color = widget->color();
           widget->lock();
           widget->setRasterOp(XorROP);
           *widget << CGAL::GREEN;
           *widget << Segment(Point(x1,y1), Point(x2,y2));
           widget->setRasterOp(old_raster);
           widget->setColor(old_color);
           widget->unlock();
	   firsttime = true;
         }
         break;
    }//endswitch
  }

  void leaveEvent(QEvent *)
  {
    if(firstpoint)
    {
      RasterOp old_raster = widget->rasterOp();//save the initial raster mode
      QColor old_color = widget->color();
      widget->lock();
        widget->setRasterOp(XorROP);
        *widget << CGAL::GREEN;
        *widget << Segment(Point(x1,y1), Point(x2,y2));
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
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      RasterOp old_raster = widget->rasterOp();//save the initial raster mode
      QColor old_color = widget->color();
      widget->setRasterOp(XorROP);
      widget->lock();
      *widget << CGAL::GREEN;
      if(!firsttime)
      *widget << Segment(Point(x1,y1),Point(x2,y2));
      *widget << Segment(Point(x1,y1),Point(x,y));
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
    oldpolicy = widget->focusPolicy();
    widget->setFocusPolicy(QWidget::StrongFocus);
    oldcursor = widget->cursor();
    widget->setCursor(cursor);
  };

  void deactivating()
  {
    widget->setFocusPolicy(oldpolicy);
    widget->setCursor(oldcursor);
    firstpoint = false;
  };

  QCursor oldcursor;
  QCursor cursor;

  FT  x1, //the X of the first point
      y1; //the Y of the first point
  FT  x2, //the old second point's X
      y2; //the old second point's Y
  bool	firstpoint, //true if the user left clicked once
        firsttime;  //true if the line is not drawn
  QWidget::FocusPolicy	oldpolicy;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_WIDGET_GET_SEGMENT_H
