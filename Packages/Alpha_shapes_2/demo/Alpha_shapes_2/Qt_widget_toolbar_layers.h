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
// file          : include/CGAL/IO/Qt_widget_toolbar_layers.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

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
