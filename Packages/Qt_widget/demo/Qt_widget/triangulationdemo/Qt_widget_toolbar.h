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
// author(s)     : Ursu Radu
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
#include <CGAL/Delaunay_triangulation_2.h>


// TODO: check if some of those includes shouldn't be in the .C file
#include <CGAL/IO/Qt_Widget.h>
#include <CGAL/IO/Qt_Widget_MovePoint.h>
#include <CGAL/IO/Qt_Widget_Get_line.h>
#include <CGAL/IO/Qt_Widget_Get_point.h>

#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qmainwindow.h>

typedef double Coord_type;
typedef CGAL::Cartesian<Coord_type>  Rp;
typedef CGAL::Delaunay_triangulation_2<Rp>  Delaunay;


namespace CGAL {

class Tools_toolbar : public QObject
{
	Q_OBJECT
public:
  Tools_toolbar(Qt_widget *w, QMainWindow *mw, Delaunay *t);
  ~Tools_toolbar()
  {
    delete maintoolbar;
  };
  QToolBar*	toolbar(){return maintoolbar;}

signals:
  void new_object(CGAL::Object);
  void was_repainted();

private slots:
  void get_new_object(CGAL::Object obj) { emit(new_object(obj)); }

  void pointtool();
  void linetool();
  void notool();
  void movepoint();

  void toggle_button();

private:
  QToolBar		*maintoolbar;
  QToolButton		*but[10];
  Qt_widget		*widget;
  int			activebutton;
  bool			is_active;
  void			setActiveButton(int i);
  void			addToolButton(QToolButton *b);
  int			nr_of_buttons;
	
  Delaunay			    *dt;
  CGAL::Qt_widget_get_line<Rp>	    linebut;
  CGAL::Qt_widget_get_point<Rp>	    pointbut;
  CGAL::Qt_widget_movepoint<Delaunay> movepointbut;
};//end class

};//end namespace

#endif
