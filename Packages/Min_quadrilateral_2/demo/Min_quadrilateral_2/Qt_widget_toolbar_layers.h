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
// file          : Qt_widget_toolbar_layers.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WIDGET_TOOLBAR_LAYERS_H
#define CGAL_QT_WIDGET_TOOLBAR_LAYERS_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Qt_widget.h>
#include <qmainwindow.h>
#include "min_quadrilateral_layers.h"

#include <qobject.h>
#include <qmainwindow.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qstatusbar.h>
#include <qbuttongroup.h>

typedef double Coord_type;
typedef CGAL::Cartesian<Coord_type>  Rp;
typedef Rp::Point_2   Point;

class Layers_toolbar : public QToolBar
{
  Q_OBJECT
public:
  Layers_toolbar(CGAL::Qt_widget *w, QMainWindow *mw, std::list<Point> *l_of_p);
  ~Layers_toolbar()
  {
    delete showP;
    delete showPL;
    delete showLS;
    delete showR;
  };
  
private:
  QToolButton     *but[10];
  CGAL::Qt_widget *widget;
  QMainWindow     *window;
  QButtonGroup    *button_group;
  int             nr_of_buttons;
	
  Qt_layer_show_points<Rp>          *showP;
  Qt_layer_show_parallelogram<Rp>   *showPL;
  Qt_layer_show_strip<Rp>           *showLS;
  Qt_layer_show_rectangle<Rp>       *showR;
};//end class

#endif
