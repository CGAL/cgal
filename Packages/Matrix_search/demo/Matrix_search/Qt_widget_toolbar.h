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
// file          : demo/Qt_widget/Max_k-gon/Qt_widget_toolbar.h
// package       : Qt_widget
// author(s)     : Ursu Radu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


#ifndef CGAL_QT_WIDGET_TOOLBAR_H
#define CGAL_QT_WIDGET_TOOLBAR_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>


// TODO: check if some of those includes shouldn't be in the .C file
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_get_point.h>
#include "Qt_widget_move_list_point.h"

#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qbuttongroup.h>
#include <qmainwindow.h>

typedef CGAL::Cartesian<double>	  Rp;
typedef Rp::Point_2		  Point;

class Tools_toolbar : public QToolBar
{
	Q_OBJECT
public:
  Tools_toolbar(CGAL::Qt_widget *w, QMainWindow *mw, std::list<Point> *l1);
  ~Tools_toolbar(){};
signals:
  void new_object(CGAL::Object);

private:
  QToolButton      *but[10];
  QButtonGroup     *button_group;
  CGAL::Qt_widget  *widget;
  int              nr_of_buttons;


  CGAL::Qt_widget_get_point<Rp>	pointbut;
  Qt_widget_move_list_point<Rp> move_deletebut;
};//end class

#endif
