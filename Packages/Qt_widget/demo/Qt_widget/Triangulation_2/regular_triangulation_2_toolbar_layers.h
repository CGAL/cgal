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
// file          : triangulation_2_toolbar_layers.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

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
