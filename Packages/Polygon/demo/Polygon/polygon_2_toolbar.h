// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
// This software and related documentation are part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation are provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind. 
// ----------------------------------------------------------------------
// file          : partition_2_toolbar.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// coordinator   : Laurent Rineau
//
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================


#ifndef CGAL_POLYGON_TOOLBAR_H
#define CGAL_POLYGON_TOOLBAR_H


#include "cgal_types.h"
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_get_simple_polygon.h>
#include <CGAL/IO/Qt_widget_get_polygon.h>
#include <CGAL/IO/Qt_widget_get_point.h>

#include <qobject.h>
#include <qmainwindow.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qbuttongroup.h>

class Polygon_toolbar : public QToolBar
{
  Q_OBJECT
public:
  Polygon_toolbar(CGAL::Qt_widget *w, QMainWindow *mw);
	
private:
  QToolButton     *but[10];
  CGAL::Qt_widget *widget;
  QButtonGroup    *button_group;
  int             nr_of_buttons;
	
  CGAL::Qt_widget_get_simple_polygon<Cgal_Polygon> getsimplepoly;
  CGAL::Qt_widget_get_polygon<Cgal_Polygon>        getpoly;
  CGAL::Qt_widget_get_point<K>                     getpoint;
};//end class

#endif
