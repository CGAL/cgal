// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
// This software and related documentation are part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation are provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind. 
// ----------------------------------------------------------------------
// file          : partition_2_toolbar_layers.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// coordinator   : Laurent Rineau
//
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================

#ifndef CGAL_PARTITION_2_TOOLBAR_LAYERS_H
#define CGAL_PARTITION_2_TOOLBAR_LAYERS_H

#include "cgal_types.h"
#include <CGAL/IO/Qt_widget.h>

#include <qobject.h>
#include <qmainwindow.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qstatusbar.h>
#include <qbuttongroup.h>


template <class T> class Qt_layer_show_polygon;
template <class T> class Qt_layer_show_greene_approx;
template <class T> class Qt_layer_show_ymonotone;
template <class T> class Qt_layer_show_optimal_convex;
template <class T> class Qt_layer_show_polygon_points;

class Layers_toolbar : public QToolBar
{
  Q_OBJECT
public:
  Layers_toolbar(CGAL::Qt_widget *w, QMainWindow *mw, Cgal_Polygon *p);
  ~Layers_toolbar();
private:
  QToolButton         *but[10];
  CGAL::Qt_widget     *widget;
  QMainWindow         *window;
  QButtonGroup        *button_group;
  int                 nr_of_buttons;

  Qt_layer_show_polygon <Cgal_Polygon>        *showP;
  Qt_layer_show_greene_approx <Cgal_Polygon > *showGA;
  Qt_layer_show_ymonotone <Cgal_Polygon>      *showYM;
  Qt_layer_show_optimal_convex <Cgal_Polygon> *showOC;
  Qt_layer_show_polygon_points <Cgal_Polygon> *showPP;

};//end class

#endif
