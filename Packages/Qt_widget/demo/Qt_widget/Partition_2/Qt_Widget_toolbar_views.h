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
// file          : include/CGAL/IO/Qt_Window_toolbar_views.h
// package       : QT_window
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WINDOW_TOOLBAR_VIEWS_H
#define CGAL_QT_WINDOW_TOOLBAR_VIEWS_H

#include <list>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Partition_traits_2.h>

#include <CGAL/IO/Qt_Widget.h>

#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qstatusbar.h>

typedef CGAL::Cartesian<double>                           K;
typedef CGAL::Partition_traits_2<K>                       Traits;
typedef Traits::Point_2                                   Point_2;
typedef Traits::Polygon_2                                 Polygon;

namespace CGAL {

class Qt_view_mouse_coordinates;
template <class T> class Qt_view_show_polygon;
template <class T> class Qt_view_show_greene_approx;
template <class T> class Qt_view_show_ymonotone;
template <class T> class Qt_view_show_optimal_convex;
template <class T> class Qt_view_show_polygon_points;

class Views_toolbar : public QObject
{
  Q_OBJECT
public:
  Views_toolbar(Qt_widget *w, QMainWindow *mw, Polygon *p);
  ~Views_toolbar()
  {
    delete maintoolbar;
  };
  QToolBar* toolbar(){return maintoolbar;};

signals:
  void new_object(CGAL::Object);
		
private slots:
  void show_coordinates();
  void show_polygon();
  void show_greene_approx();
  void show_ymonotone();
  void show_optimal();
  void show_points();

private:
  QToolBar		*maintoolbar;
  QToolButton		*but[10];
  Qt_widget		*widget;
  QMainWindow		*window;	
  int			nr_of_buttons;

  CGAL::Qt_view_mouse_coordinates		*showMC;
  CGAL::Qt_view_show_polygon <Polygon>		*showP;
  CGAL::Qt_view_show_greene_approx <Polygon >	*showGA;
  CGAL::Qt_view_show_ymonotone <Polygon>	*showYM;
  CGAL::Qt_view_show_optimal_convex <Polygon>	*showOC;
  CGAL::Qt_view_show_polygon_points <Polygon>	*showPP;

};//end class

};//end namespace

#endif
