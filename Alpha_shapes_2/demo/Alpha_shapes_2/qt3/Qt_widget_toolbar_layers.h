// Copyright (c) 2002  INRIA Sophia-Antipolis (France).
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


//Qt_widget
#include "cgal_types.h"
#include <CGAL/IO/Qt_widget.h>

//Qt_widget_layer
#include "alpha_shapes_2_layers.h"

//Qt
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
  Layers_toolbar(CGAL::Qt_widget *w, QMainWindow *mw,
                 Delaunay *t, Alpha_shape *a, QImage *i);

private:
  QToolButton     *but[10];
  QButtonGroup    *button_group;
  CGAL::Qt_widget *widget;
  QMainWindow     *window;
  int             nr_of_buttons;

  Qt_layer_show_triangulation < Delaunay >  *showT;
  Qt_layer_show_voronoi < Delaunay >        *showV;
  Qt_layer_show_points < Delaunay >         *showP;
  Qt_layer_show_alpha_shape                 *showA;
  Qt_layer_show_image                       *showI;
};//end class

#endif
