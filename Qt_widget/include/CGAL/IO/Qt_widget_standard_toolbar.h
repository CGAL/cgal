// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Radu Ursu

#ifndef CGAL_QT_WIDGET_STANDARD_TOOLBAR_H
#define CGAL_QT_WIDGET_STANDARD_TOOLBAR_H

#include <CGAL/basic.h>

#include <qtoolbar.h>
#include <qbuttongroup.h>
#include <qtoolbutton.h>

namespace CGAL {

class Qt_widget;
class Qt_widget_history;

class Qt_widget_standard_toolbar : public QToolBar
{
	Q_OBJECT
public:
  Qt_widget_standard_toolbar(Qt_widget *w,
			     QMainWindow *parent = 0,
			     const char* name = 0);

  Qt_widget_standard_toolbar(Qt_widget *w,
			     QMainWindow *mw,
			     QWidget* parent,
			     bool newLine = true,
			     const char* name = 0);

  ~Qt_widget_standard_toolbar() { delete button_group; }

  // CGAL-2.4 compatibility
  QToolBar*	toolbar(){return this;}

public slots:
  void back();
  void forward();
  void clear_history();

private slots:
  void zoomin();
  void zoomout();
  void group_clicked(int i);

private:
  void fill_toolbar(QMainWindow *mw);

private:
  Qt_widget          *widget;
  Qt_widget_history  *history;
  QButtonGroup*      button_group;
  // this group has no parent and is destroyed manually in the
  // destructor

  QToolButton* nolayerBt;
};//end class

} //end namespace

#endif //CGAL_QT_WIDGET_STANDARD_TOOLBAR_H
