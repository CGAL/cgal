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
// file          : include/CGAL/IO/Qt_widget_toolbar_layers.h
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

#include <list>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Partition_traits_2.h>

#include <CGAL/IO/Qt_widget.h>

#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qstatusbar.h>
#include <qbuttongroup.h>

typedef CGAL::Cartesian<CGAL::MP_Float>                   K;
typedef CGAL::Partition_traits_2<K>                       Traits;
typedef Traits::Point_2                                   Point_2;
typedef Traits::Polygon_2                                 Polygon;

namespace CGAL {

class Qt_layer_mouse_coordinates;
template <class T> class Qt_layer_show_polygon;
template <class T> class Qt_layer_show_greene_approx;
template <class T> class Qt_layer_show_ymonotone;
template <class T> class Qt_layer_show_optimal_convex;
template <class T> class Qt_layer_show_polygon_points;

class Layers_toolbar : public QObject
{
  Q_OBJECT
public:
  Layers_toolbar(Qt_widget *w, QMainWindow *mw, Polygon *p);
  ~Layers_toolbar();
  QToolBar* toolbar(){return maintoolbar;};

signals:
  void new_object(CGAL::Object);
		
private slots:
  void redraw_win(int);
private:
  QToolBar		  *maintoolbar;
  QToolButton		*but[10];
  Qt_widget		  *widget;
  QMainWindow		*window;
  QButtonGroup  *button_group;
  int			nr_of_buttons;

  CGAL::Qt_layer_mouse_coordinates		          *showMC;
  CGAL::Qt_layer_show_polygon <Polygon>		      *showP;
  CGAL::Qt_layer_show_greene_approx <Polygon >	*showGA;
  CGAL::Qt_layer_show_ymonotone <Polygon>	      *showYM;
  CGAL::Qt_layer_show_optimal_convex <Polygon>	*showOC;
  CGAL::Qt_layer_show_polygon_points <Polygon>	*showPP;

};//end class

};//end namespace

#endif
