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
// file          : include/CGAL/IO/Qt_Window_toolbar.h
// package       : QT_window
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


#ifndef CGAL_QT_WINDOW_TOOLBAR_H
#define CGAL_QT_WINDOW_TOOLBAR_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Partition_traits_2.h>
#include <list>


// TODO: check if some of those includes shouldn't be in the .C file
#include <CGAL/IO/Qt_Widget.h>
#include <CGAL/IO/Qt_Widget_Zoom.h>
#include <CGAL/IO/Qt_Widget_Zoomrect.h>
#include <CGAL/IO/Qt_Widget_Handtool.h>
#include <CGAL/IO/Qt_Widget_Get_simple_polygon.h>


#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qmainwindow.h>

typedef CGAL::Cartesian<double>                           Rp;
typedef CGAL::Partition_traits_2<Rp>			  Traits;
typedef Traits::Point_2                                   Point_2;
typedef Traits::Polygon_2                                 Polygon;


namespace CGAL {

class Tools_toolbar : public QObject
{
  Q_OBJECT
public:
  Tools_toolbar(Qt_widget *w, QMainWindow *mw);
  ~Tools_toolbar()
  {
    delete maintoolbar;
  };
  QToolBar*	toolbar(){return maintoolbar;}

signals:
  void new_object(CGAL::Object);
  void was_repainted();

private slots:
  void toolregion();
  void zoomin();
  void zoomout();
  void zoominrect();
  void notool();
  void handtool();
  void get_s_polygon();
	
	
private:
  QToolBar				*maintoolbar;
  QToolButton				*but[10];
  Qt_widget				*widget;
  int					activebutton;
  bool					is_active;
  void					setActiveButton(int i);
  void					addToolButton(QToolButton *b);
  int					nr_of_buttons;
	
  CGAL::Qt_widget_zoom			zoombut;
  CGAL::Qt_widget_zoomrect		zoomrectbut;
  CGAL::Qt_widget_handtool<Rp>		handtoolbut;
  CGAL::Qt_widget_get_simple_polygon<Polygon> getsimplebut;
};//end class

};//end namespace

#endif
