// Copyright (c) 1997-2002  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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


#ifndef CGAL_QT_WIDGET_TOOLBAR_LAYERS_H
#define CGAL_QT_WIDGET_TOOLBAR_LAYERS_H

#include "constrained_cgal_types.h"
#include <CGAL/IO/Qt_widget.h>
#include "constrained_delaunay_triangulation_2_layers.h"


#include <qobject.h>
#include <qmainwindow.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qstatusbar.h>
#include <qbuttongroup.h>


class Layers_toolbar : public QToolBar
{
  Q_OBJECT
public:
  Layers_toolbar(CGAL::Qt_widget *w, QMainWindow *mw, CDT *t);
  ~Layers_toolbar()
  {
    delete showT;
    delete showP;
    delete showC;
  };
private:
  QToolButton      *but[10];
  CGAL::Qt_widget  *widget;
  QMainWindow      *window;
  QButtonGroup     *button_group;
  int              nr_of_buttons;

  Qt_layer_show_triangulation < CDT>  *showT;
  Qt_layer_show_points < CDT >        *showP;
  Qt_layer_show_constraints < CDT >   *showC;
};//end class

#endif
