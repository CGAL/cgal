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
// file          : include/CGAL/IO/Qt_widget_show_mouse_coordinates.h
// package       : Qt_widget
// author(s)     : Ursu Radu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WIDGET_SHOW_MOUSE_COORDINATES_H
#define CGAL_QT_WIDGET_SHOW_MOUSE_COORDINATES_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <qobject.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qstring.h>

namespace CGAL {

class Qt_widget_show_mouse_coordinates : public Qt_widget_layer
{
public:
	
  Qt_widget_show_mouse_coordinates(QMainWindow &mw) : qmw(mw){};
  void draw(){};
  void mouseMoveEvent(QMouseEvent *e)
  {
    QString s("x=%1 y=%2");
    double xcoord, ycoord;
    widget->x_real(e->x(), xcoord);
    widget->y_real(e->y(), ycoord);
    qmw.statusBar()->message(s.arg(xcoord, -20, 'g', 15).
			     arg(ycoord, -20,'g', 15));

  };
private:
  void deactivating(){
    qmw.statusBar()->clear();
  }
  QMainWindow	&qmw;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_WIDGET_GET_SEGMENT_H
