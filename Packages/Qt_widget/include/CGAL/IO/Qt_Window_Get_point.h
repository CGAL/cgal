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
// file          : include/CGAL/IO/Qt_Window_Get_point.h
// package       : QT_window
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WINDOW_GET_POINT_H
#define CGAL_QT_WINDOW_GET_POINT_H

#include <CGAL/IO/Qt_Window.h>
#include <CGAL/IO/Qt_Window_tool.h>

#ifndef CGAL_QT_WINDOW_GET_POINT_BUTTON
#define CGAL_QT_WINDOW_GET_POINT_BUTTON Qt::LeftButton
#endif

namespace CGAL {

template <class R>
class Qt_widget_get_point : public Qt_widget_tool
{
  //  Q_OBJECT
public:
  typedef Point_2<R> Point;
  typedef typename R::FT FT;

  Qt_widget_get_point() {};

private:
  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == CGAL_QT_WINDOW_GET_POINT_BUTTON)
      {
	FT
	  x=static_cast<FT>(widget->x_real(e->x())),
	  y=static_cast<FT>(widget->y_real(e->y()));
	emit new_object(make_object(Point(x,y)));
      }
  };
};

} // namespace CGAL

#endif // CGAL_QT_WINDOW_GET_POINT_H
