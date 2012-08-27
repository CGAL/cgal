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
// Author(s)     : Laurent Rineau

#ifndef CGAL_QT_WIDGET_GET_POINT_H
#define CGAL_QT_WIDGET_GET_POINT_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <qcursor.h>

#ifndef CGAL_QT_WIDGET_GET_POINT_BUTTON
#define CGAL_QT_WIDGET_GET_POINT_BUTTON Qt::LeftButton
#endif

namespace CGAL {

template <class R>
class Qt_widget_get_point : public Qt_widget_layer
{
public:
  typedef typename R::Point_2	Point;
  typedef typename R::FT	FT;
  
  Qt_widget_get_point(const QCursor c=QCursor(Qt::crossCursor),
		      QObject* parent = 0, const char* name = 0) :
    Qt_widget_layer(parent, name), cursor(c) {};
  
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
       && is_pure(e->state()))
    {
      FT x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      widget->new_object(make_object(Point(x, y)));
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
  };

  QCursor cursor;
  QCursor oldcursor;
};

} // namespace CGAL

#endif // CGAL_QT_WIDGET_GET_POINT_H
