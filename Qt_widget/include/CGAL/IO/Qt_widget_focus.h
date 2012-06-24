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

#ifndef CGAL_QT_WIDGET_FOCUS_H
#define CGAL_QT_WIDGET_FOCUS_H

#include <CGAL/IO/pixmaps/focus1.xpm>
#include <CGAL/IO/pixmaps/focus1_mask.xpm>
#include <CGAL/IO/pixmaps/focus2.xpm>
#include <CGAL/IO/pixmaps/focus2_mask.xpm>
#include <CGAL/IO/pixmaps/focus3.xpm>
#include <CGAL/IO/pixmaps/focus3_mask.xpm>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <qcolor.h>
#include <qtimer.h>
#include <qpixmap.h>
#include <qbitmap.h>
#include <qcursor.h>


namespace CGAL {

class Qt_widget_focus : public Qt_widget_layer
{
private:
  QPixmap *mouse_ico1,
          *mouse_ico2,
          *mouse_ico3;
  QCursor *cursor1,
          *cursor2,
          *cursor3;
  QBitmap cb, cm;
  int	  cycle;
  QCursor oldcursor;

public:
  Qt_widget_focus(QObject* parent = 0, const char* name = 0) 
    : Qt_widget_layer(parent, name), cycle(0)
  {
	mouse_ico1 = new QPixmap( (const char**)focus1_xpm);
	mouse_ico2 = new QPixmap( (const char**)focus2_xpm);
	mouse_ico3 = new QPixmap( (const char**)focus3_xpm);
	
	QPixmap *mouse_ico_mask1 = new QPixmap((const char**)focus1_mask_xpm);	
	QPixmap *mouse_ico_mask2 = new QPixmap((const char**)focus2_mask_xpm);
	QPixmap *mouse_ico_mask3 = new QPixmap((const char**)focus3_mask_xpm);

	cb = *mouse_ico1; cm = *mouse_ico_mask1;
	mouse_ico1->setMask(cm);
	cursor1 = new QCursor(*mouse_ico1);

	cb = *mouse_ico2; cm = *mouse_ico_mask2;
	mouse_ico2->setMask(cm);
	cursor2 = new QCursor(*mouse_ico2);

	cb = *mouse_ico3; cm = *mouse_ico_mask3;
	mouse_ico3->setMask(cm);
	cursor3 = new QCursor(*mouse_ico3);
  };

  void timerEvent( QTimerEvent *)
  {
	switch(cycle){
	case 1:
		widget->setCursor(*cursor1);
		cycle++;
		break;
	case 2:
		widget->setCursor(*cursor2);
		cycle++;
		break;
	case 3:
		widget->setCursor(*cursor3);
		cycle=1;
		break;
	}
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
    if(e->button() == Qt::LeftButton
       && is_pure(e->state()))
    {
      double x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);	
      widget->set_center(x, y);
    }
  };
  void deactivating()
  {    
    widget->setCursor(oldcursor);
	  killTimers();
  };

  void activating()
  {
    oldcursor = widget->cursor();
	  startTimer( 200 );
	  cycle = 1;

  };
};//end class 

} // namespace CGAL

#endif // CGAL_QT_WIDGET_FOCUS_H
