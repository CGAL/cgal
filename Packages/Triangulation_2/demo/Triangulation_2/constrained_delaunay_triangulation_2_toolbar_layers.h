// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
// This software and related documentation are part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation are provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind. 
//
// Every use of CGAL requires a license. 
//
// Academic research and teaching license
// - For academic research and teaching purposes, permission to use and copy
//   the software and its documentation is hereby granted free of charge,
//   provided that it is not a component of a commercial product, and this
//   notice appears in all copies of the software and related documentation. 
//
// Commercial licenses
// - Please check the CGAL web site http://www.cgal.org/index2.html for 
//   availability.
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// file          : constrained_delaunay_triangulation_2_toolbar_layers.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau
//
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================

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
