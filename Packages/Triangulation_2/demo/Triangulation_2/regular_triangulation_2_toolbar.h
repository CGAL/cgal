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
// file          : include/CGAL/IO/Qt_widget_toolbar.h
// package       : Qt_widget
// author(s)     : Ursu Radu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_REGULAR_TRIANGULATION_2_TOOLBAR_H
#define CGAL_REGULAR_TRIANGULATION_2_TOOLBAR_H

#include "regular_cgal_types.h"
#include <CGAL/IO/Qt_widget.h>
#include "triangulation_2_edit_vertex.h"
#include <CGAL/IO/Qt_widget_get_circle.h>
#include <CGAL/IO/Qt_widget_get_point.h>

#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qmainwindow.h>
#include <qbuttongroup.h>


class Tools_toolbar : public QToolBar
{
  Q_OBJECT
public:
  Tools_toolbar(CGAL::Qt_widget *w, QMainWindow *mw, Regular_triangulation *t);
private:
  QToolButton        *but[10];
  CGAL::Qt_widget    *widget;
  QButtonGroup       *button_group;
  int                activebutton;
  bool               is_active;  
  int                nr_of_buttons;
	
  CGAL::Qt_widget_get_circle<Rp>       input_circle_layer;
  CGAL::Qt_widget_get_point<Rp>        input_point_layer;
  triangulation_2_edit_weightedpoint<Regular_triangulation>
                                       edit_vertex_layer;
};//end class


#endif
