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
// file          : include/CGAL/IO/Qt_widget_standard_toolbar.h
// package       : Qt_widget
// author(s)     : Radu Ursu
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
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_focus.h>
#include <CGAL/IO/Qt_widget_zoomrect.h>
#include <CGAL/IO/Qt_widget_handtool.h>


#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qmainwindow.h>
#include <qbuttongroup.h>

namespace CGAL {

class Qt_widget_standard_toolbar : public QObject
{
	Q_OBJECT
public:
  Qt_widget_standard_toolbar(Qt_widget *w, QMainWindow *mw);
  QToolBar*	toolbar(){return maintoolbar;}


private slots:
  void zoomin();
  void zoomout();
  void back();
  void forward();
  
private:
  QToolBar		  *maintoolbar;
  QToolButton	  *but[10];
  Qt_widget		  *widget;
  QButtonGroup  *button_group;
  int			      nr_of_buttons;
	
  CGAL::Qt_widget_focus         zoombut;
  CGAL::Qt_widget_zoomrect	    zoomrectbut;
  CGAL::Qt_widget_handtool	    handtoolbut;
};//end class

};//end namespace

#endif //CGAL_QT_WIDGET_STANDARD_TOOLBAR_H
