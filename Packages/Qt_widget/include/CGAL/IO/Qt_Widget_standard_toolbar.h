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


#ifndef CGAL_QT_WIDGET_STANDARD_TOOLBAR_H
#define CGAL_QT_WIDGET_STANDARD_TOOLBAR_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>


// TODO: check if some of those includes shouldn't be in the .C file
#include <CGAL/IO/Qt_Widget.h>
#include <CGAL/IO/Qt_Widget_Zoom.h>
#include <CGAL/IO/Qt_Widget_Zoomrect.h>
#include <CGAL/IO/Qt_Widget_Handtool.h>


#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qmainwindow.h>

typedef double Coord_type;
typedef CGAL::Cartesian<Coord_type>  Rp;

namespace CGAL {

class Standard_toolbar : public QObject
{
	Q_OBJECT
public:
  Standard_toolbar(Qt_widget *w, QMainWindow *mw);
  ~Standard_toolbar()
  {
    delete maintoolbar;
  };
  QToolBar*	toolbar(){return maintoolbar;}

signals:
  void new_object(CGAL::Object);
  void was_repainted();

private slots:
  void get_new_object(CGAL::Object obj) { emit(new_object(obj)); }

  void toggle_button();

  void toolregion();
  void zoomin();
  void zoomout();
  void zoominrect();
  void notool();
  void handtool();
	
private:
  QToolBar		*maintoolbar;
  QToolButton		*but[10];
  Qt_widget		*widget;
  int			activebutton;
  bool			is_active;
  void			setActiveButton(int i);
  int			nr_of_buttons;
	
  CGAL::Qt_widget_zoom		    zoombut;
  CGAL::Qt_widget_zoomrect	    zoomrectbut;
  CGAL::Qt_widget_handtool<Rp>	    handtoolbut;
};//end class

};//end namespace

#endif
