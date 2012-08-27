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

#ifndef CGAL_QT_WIDGET_GET_ISO_RECTANGLE_H
#define CGAL_QT_WIDGET_GET_ISO_RECTANGLE_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <qrect.h>
#include <qcursor.h>

#ifndef CGAL_QT_LEFT_BUTTON
#define CGAL_QT_LEFT_BUTTON Qt::LeftButton
#endif


namespace CGAL {
template <class T>
class Qt_widget_get_iso_rectangle : public Qt_widget_layer
{
private:
  QCursor cursor;
  QCursor oldcursor;

public:
  int                                 first_x, first_y, x2, y2;
  bool                                widgetrepainted;
  bool                                on_first;
  QWidget::FocusPolicy                oldpolicy;
  typedef typename T::Iso_rectangle_2 Iso_rectangle_2;
  typedef typename T::RT              RT;

  Qt_widget_get_iso_rectangle(const QCursor
			      c=QCursor(Qt::crossCursor),
			      QObject* parent = 0, const char* name = 0)
    : Qt_widget_layer(parent, name),  cursor(c), widgetrepainted(true),
      on_first(false) {};
  void draw(){
    widgetrepainted = true;
  };

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
    if(e->button() == CGAL_QT_LEFT_BUTTON
       && is_pure(e->state()))
    {
      if (!on_first)
      {
        first_x = e->x();
        first_y = e->y();
        on_first = true;
      } else {
          if((e->x() != first_x) && (e->y() != first_y)) {
          RT x, y, xfirst2, yfirst2;
          widget->x_real(e->x(), x);
          widget->y_real(e->y(), y);
          widget->x_real(first_x, xfirst2);
          widget->y_real(first_y, yfirst2);
          RT xmin, xmax, ymin, ymax;
          if(x < xfirst2) {xmin = x; xmax = xfirst2;}
          else {xmin = xfirst2; xmax = x;};
          if(y < yfirst2) {ymin = y; ymax = yfirst2;}
          else {ymin = yfirst2; ymax = y;};

          widget->new_object(
                     make_object(Iso_rectangle_2(xmin, ymin, 
                                                 xmax, ymax)));
          on_first = false;
          widgetrepainted = true;
        }
      }
    }
  };


  void mouseMoveEvent(QMouseEvent *e)
  {
    if(on_first)
    {
      int x = e->x();
      int y = e->y();
      *widget << noFill;
      RasterOp old = widget->rasterOp();	//save the initial raster mode
      QColor old_color=widget->color();
      widget->setRasterOp(XorROP);
      widget->lock();
      widget->setColor(Qt::green);
      if(!widgetrepainted)
        widget->get_painter().drawRect(first_x, first_y, 
                                       x2 - first_x, y2 - first_y);
      widget->get_painter().drawRect(first_x, first_y, x - first_x,
                                     y - first_y);
      widget->unlock();
      widget->setColor(old_color);
      widget->setRasterOp(old);

      //save the last coordinates to redraw the screen
      x2 = x;
      y2 = y;
      widgetrepainted = false;
    }
  };

  void keyPressEvent(QKeyEvent *e)
  {
    switch ( e->key() ) {
      case Key_Escape:			// key_escape
         if (on_first)
         {
           widget->lock();
           *widget << noFill;
           RasterOp old = widget->rasterOp();	//save the initial raster mode
           QColor old_color=widget->color();
           widget->setRasterOp(XorROP);
           *widget << CGAL::GREEN;
           if(!widgetrepainted)
             widget->get_painter().drawRect(first_x, first_y, 
                                       x2 - first_x, y2 - first_y);
           widget->setColor(old_color);
           widget->setRasterOp(old);
           widgetrepainted = true;

           widget->unlock();
	   on_first = false;
         }
         break;
    }//endswitch
  }

  void leaveEvent(QEvent *)
  {
    if (on_first)
    {
      widget->lock();
      *widget << noFill;
      RasterOp old = widget->rasterOp();	//save the initial raster mode
      QColor old_color=widget->color();
      widget->setRasterOp(XorROP);
      *widget << CGAL::GREEN;
      if(!widgetrepainted)
        widget->get_painter().drawRect(first_x, first_y, 
                                       x2 - first_x, y2 - first_y);
      widget->setColor(old_color);
      widget->setRasterOp(old);
      widgetrepainted = true;

      widget->unlock();
    }
  }

  void activating()
  {
    oldpolicy = widget->focusPolicy();
    widget->setFocusPolicy(QWidget::StrongFocus);
    oldcursor = widget->cursor();
    widget->setCursor(cursor);
    widgetrepainted = true;
  };

  void deactivating()
  {
    widget->setCursor(oldcursor);
    widget->setFocusPolicy(oldpolicy);
    on_first = false;
  };
};//end class 

} // namespace CGAL


#endif
