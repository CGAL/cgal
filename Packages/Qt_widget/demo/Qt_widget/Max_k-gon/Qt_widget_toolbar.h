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


// TODO: check if some of those includes shouldn't be in the .C file
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_get_point.h>

#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qmainwindow.h>

typedef CGAL::Cartesian<double>	  Rp;
typedef Rp::Point_2		  Point;

namespace CGAL {

class Tools_toolbar : public QObject
{
	Q_OBJECT
public:
  Tools_toolbar(Qt_widget *w, QMainWindow *mw, std::list<Point> &l1);
  ~Tools_toolbar()
  {
    delete maintoolbar;
  };
  QToolBar*	toolbar(){return maintoolbar;}

signals:
  void new_object(CGAL::Object);


private slots:
  void get_new_object(CGAL::Object obj) { emit(new_object(obj)); }

  void pointtool();
  void notool();
  
  void toggle_button();

private:
  QToolBar		*maintoolbar;
  QToolButton		*but[10];
  Qt_widget		*widget;
  int			activebutton;
  bool			is_active;
  void			setActiveButton(int i);
  int			nr_of_buttons;


  CGAL::Qt_widget_get_point<Rp>	    pointbut;
};//end class

};//end namespace

#endif
