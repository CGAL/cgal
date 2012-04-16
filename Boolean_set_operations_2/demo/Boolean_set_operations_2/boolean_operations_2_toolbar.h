// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef BOOLEAN_OPERATIONS_2_TOOLBAR_H
#define BOOLEAN_OPERATIONS_2_TOOLBAR_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>



#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_get_point.h>
//#include <CGAL/IO/Qt_widget_get_simple_polygon.h>
#include "Qt_widget_get_circ_polygon.h"

#include <CGAL/IO/Qt_widget_get_circle.h>
#include <CGAL/IO/pixmaps/circle.xpm>


#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qbuttongroup.h>
#include <qmainwindow.h>
#include <qcursor.h>
#include <qradiobutton.h>
#include <qvbuttongroup.h>



#include "typedefs.h"
#include "Qt_widget_locate_layer.h"

class MyWindow;

class Tools_toolbar : public QToolBar
{
  Q_OBJECT
public:
	//mw type was QMainWindow 
	Tools_toolbar(CGAL::Qt_widget *w, MyWindow *mw);
  ~Tools_toolbar(){};

  void deactivate()
  {
    getsimplebut.deactivate();
    getcirclebut.deactivate();
    but[0]->toggle(); // toggle the 'deactivate layer button'
  }

  void reset()
  {
    locatebut->reset();
  }

private:
  QToolButton     *but[10];
  QButtonGroup    *button_group;
  CGAL::Qt_widget *widget;

  CGAL::Qt_widget_get_circ_polygon<Kernel>     getsimplebut;
  CGAL::Qt_widget_get_circle<Kernel>           getcirclebut;
  Qt_widget_locate_layer*                       locatebut;

};//end class


#endif
