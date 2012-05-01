// Copyright (c) 1997-2002  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Radu Ursu

#ifndef CGAL_TRIANGULATION_2_TOOLBAR_LAYERS_H
#define CGAL_TRIANGULATION_2_TOOLBAR_LAYERS_H

#include "regular_cgal_types.h"
#include <CGAL/IO/Qt_widget.h>
#include "regular_triangulation_2_layers.h"

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
  Layers_toolbar(CGAL::Qt_widget *w, QMainWindow *mw, Regular_triangulation *t);
  ~Layers_toolbar()
  {
    delete showT;
    delete showV;
    delete showP;
  };
private:
  QToolButton     *but[10];
  CGAL::Qt_widget *widget;
  QMainWindow     *window;
  Regular_triangulation
                  *dt;
  QButtonGroup    *button_group;
  int             nr_of_buttons;

  Qt_layer_show_triangulation < Regular_triangulation > *showT;
  Qt_layer_show_voronoi < Regular_triangulation >       *showV;
  Qt_layer_show_points < Regular_triangulation >        *showP;
};//end class

#endif
