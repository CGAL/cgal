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

#ifndef CGAL_QT_WIDGET_HANDTOOL_H
#define CGAL_QT_WIDGET_HANDTOOL_H

#include <cstdio>
#include <CGAL/IO/pixmaps/hand.xpm>
#include <CGAL/IO/pixmaps/holddown.xpm>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <qrect.h>
#include <qcursor.h>

#ifndef CGAL_QT_WIDGET_GET_POINT_BUTTON
#define CGAL_QT_WIDGET_GET_POINT_BUTTON Qt::LeftButton
#endif


namespace CGAL {

class Qt_widget_handtool : public Qt_widget_layer
{
public:
  Qt_widget_handtool(QObject* parent = 0, const char* name = 0)
    : Qt_widget_layer(parent, name), wasrepainted(TRUE), on_first(FALSE){};

private:
  QCursor oldcursor;

  void draw(){
    wasrepainted = TRUE;
  };

  void timerEvent( QTimerEvent *)
  {
    if(on_first)
      widget->setCursor(QCursor( 
              QPixmap( (const char**)holddown_xpm)));
    else
      widget->setCursor(QCursor( 
              QPixmap( (const char**)hand_xpm)));
  }

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
       && is_pure(e->state()))
    {
      widget->setCursor(QCursor( QPixmap( (const char**)holddown_xpm)));
      if (!on_first){
	      first_x = e->x();
	      first_y = e->y();
	      on_first = TRUE;
      }	
    }
  };

  void mouseReleaseEvent(QMouseEvent *e)
  {
    if(e->button() == CGAL_QT_WIDGET_GET_POINT_BUTTON
       && is_pure(e->state()))
    {
      widget->setCursor(QCursor( QPixmap( (const char**)hand_xpm)));
      double x, y, xfirst2, yfirst2;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      widget->x_real(first_x, xfirst2);
      widget->y_real(first_y, yfirst2);
			
      double distx, disty;
      distx = xfirst2 - x;
      disty = yfirst2 - y;
      widget->move_center(distx, disty);
      on_first = FALSE;
    }
  }
  void mouseMoveEvent(QMouseEvent *e)
  {
    char tempc1[130], tempc2[40];
    double xcoord, ycoord;
    if(on_first)
    {
      int x = e->x();
      int y = e->y();
      //save the initial raster mode
      RasterOp old = widget->rasterOp();	
      widget->setRasterOp(XorROP);
      widget->lock();
        *widget << CGAL::GRAY;
      if(!wasrepainted) {
        widget->x_real(x2 - first_x, xcoord);
        widget->x_real(y2 - first_y, ycoord);
        std::sprintf(tempc1, " dx=%20.6f", xcoord);
        std::sprintf(tempc2, ", dy=%20.6f", ycoord);
        strcat(tempc1, tempc2);
        widget->get_painter().drawLine(first_x, first_y, x2, y2);
        *widget << CGAL::GREEN;
        widget->get_painter().drawText(x2, y2, tempc1, 49);
        *widget << CGAL::GRAY;
      }
      widget->x_real(x - first_x, xcoord);
      widget->x_real(y - first_y, ycoord);
      std::sprintf(tempc1, " dx=%20.6f", xcoord);
      std::sprintf(tempc2, ", dy=%20.6f", ycoord);
      strcat(tempc1, tempc2);
      widget->get_painter().drawLine(first_x, first_y, x, y);
      *widget << CGAL::GREEN;
      widget->get_painter().drawText(x, y, tempc1, 49);
      widget->unlock();
      widget->setRasterOp(old);

      //save the last coordinates to redraw the screen
      x2 = x;
      y2 = y;
      wasrepainted = FALSE;
    }
  };

  void activating()
  {
    oldcursor = widget->cursor();
    widget->setCursor(QCursor( QPixmap( (const char**)hand_xpm)));
    wasrepainted = TRUE;
	  startTimer( 100 );
  };

  void deactivating()
  {
    widget->setCursor(oldcursor);
    killTimers();
  };

  int   first_x, first_y;
  int   x2, y2;
  bool	wasrepainted;
  bool	on_first;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_WIDGET_HANDTOOL_H
