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

#include <qtoolbar.h>

namespace CGAL {

class Qt_widget;

class Qt_widget_standard_toolbar : public QToolBar
{
	Q_OBJECT
public:
  Qt_widget_standard_toolbar(Qt_widget *w,
			     QMainWindow *mw,
#if QT_VERSION < 300
			     // for Qt 2.3 and before
			     QMainWindow::ToolBarDock = QMainWindow::Top,
#else
			     // from Qt 3.0
			     Dock = DockTop,
#endif
			     bool newLine = true,
			     const char* name = 0);

  QToolBar*	toolbar(){return this;} // for compatibility

private slots:
  void zoomin();
  void zoomout();
  void back();
  void forward();
  
private:
  Qt_widget	  *widget;
};//end class

};//end namespace

#endif //CGAL_QT_WIDGET_STANDARD_TOOLBAR_H
