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
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// file          : src/Qt_widget_toolbar_layers.C
// package       : Qt_widget
// author(s)     : Ursu Radu
// release       : CGAL-2.4
// release_date  : 2002, July 8
//
// coordinator   : Laurent Rineau
//
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================


#if defined CGAL_USE_QT && defined CGAL_USE_GMP

#include "Qt_widget_toolbar_layers.h"


Layers_toolbar::Layers_toolbar( CGAL::Qt_widget *w, QMainWindow *mw,
                                  Nef_polyhedron *n1, Nef_polyhedron *n2) : 
    widget(w), window(mw), nr_of_buttons(0)
  {
    showNG  = new Qt_layer_nef_gray<Nef_polyhedron>(*n1);
    showNR  = new Qt_layer_nef_red<Nef_polyhedron>(*n2);

    //set the widget
    window->statusBar();

    widget->attach(showNG);    
    widget->attach(showNR);
    
    
    maintoolbar = new QToolBar("tools", mw, QMainWindow::Top, TRUE, "Tools");

  }//end constructor

#include "Qt_widget_toolbar_layers.moc"

#endif
