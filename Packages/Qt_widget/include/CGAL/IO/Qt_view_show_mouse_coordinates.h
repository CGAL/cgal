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
// file          : include/CGAL/IO/Qt_view_show_mouse_coordinates.h
// package       : QT_window
// author(s)     : Ursu Radu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_VIEW_MOUSE_COORDINATES_H
#define CGAL_QT_VIEW_MOUSE_COORDINATES_H

#include <CGAL/IO/Qt_widget_view.h>
#include <qobject.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <cstdio>

namespace CGAL {

class Qt_view_mouse_coordinates : public Qt_widget_view
{
public:
	
  Qt_view_mouse_coordinates(QMainWindow &mw) : qmw(mw){};
  void draw_view(Qt_widget &widget){};
  void mouseMoveEvent(QMouseEvent *e, Qt_widget &widget)
  {
    char xsir[40], ysir[40], final[80];
    CGAL_CLIB_STD::sprintf(xsir, "%.15f", widget.x_real(e->x()));
    CGAL_CLIB_STD::sprintf(ysir, "%.15f", widget.y_real(e->y()));
    CGAL_CLIB_STD::sprintf(final, "x=%s  y=%s", xsir, ysir);
    qmw.statusBar()->message(final);
  };
private:
  QMainWindow	&qmw;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_WINDOW_GET_SEGMENT_H
