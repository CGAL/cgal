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
// file          : delaunay_triangulation_2_toolbar.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_DELAUNAY_TRIANGULATION_2_TOOLBAR_H
#define CGAL_DELAUNAY_TRIANGULATION_2_TOOLBAR_H

#include "cgal_types.h"
#include <CGAL/IO/Qt_widget.h>
#include "triangulation_2_edit_vertex.h"
#include <CGAL/IO/Qt_widget_get_line.h>
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
  Tools_toolbar(CGAL::Qt_widget *w, QMainWindow *mw, Delaunay *t);
  void set_line_enabled(bool t){but[2]->setEnabled(t);}
  void set_move_enabled(bool t){but[3]->setEnabled(t);}
private:
  QToolButton        *but[10];
  CGAL::Qt_widget    *widget;
  QButtonGroup       *button_group;
  int                activebutton;
  bool               is_active;
  int                nr_of_buttons;
  Delaunay           *dt;
private slots:
  void               triangulation_changed(){emit(changed());}
signals:
  void               changed();
private:	
  CGAL::Qt_widget_get_line<Rep>         input_line_layer;
  CGAL::Qt_widget_get_point<Rep>        input_point_layer;
  triangulation_2_edit_vertex<Delaunay> edit_vertex_layer;
};//end class


#endif
