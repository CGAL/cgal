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
#include <qbuttongroup.h>

namespace CGAL {

class Qt_widget;
class Qt_widget_history;

class Qt_widget_standard_toolbar : public QToolBar
{
	Q_OBJECT
public:
  Qt_widget_standard_toolbar(Qt_widget *w,
			     QMainWindow *mw = 0,
			     const char* name = 0);

  Qt_widget_standard_toolbar(Qt_widget *w,
			     QMainWindow *mw,
			     QWidget* parent,
			     bool newLine = true,
			     const char* name = 0);

  ~Qt_widget_standard_toolbar() { delete button_group; }

  // CGAL-2.4 compatibility
  QToolBar*	toolbar(){return this;}

private slots:
  void zoomin();
  void zoomout();
  void back();
  void forward();
private:
  void fill_toolbar(QMainWindow *mw);
  
private:
  Qt_widget	  *widget;
  Qt_widget_history* history;
  QButtonGroup* button_group; 
    // this button has no parent and is destroyed manually in the
    // destructor

};//end class

};//end namespace

#endif //CGAL_QT_WIDGET_STANDARD_TOOLBAR_H
