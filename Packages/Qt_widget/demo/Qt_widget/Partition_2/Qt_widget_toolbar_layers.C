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
// release_date  : 2002, May 16
//
// coordinator   : Laurent Rineau
//
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================


#ifdef CGAL_USE_QT

#include "Qt_widget_toolbar_layers.h"

#include "Qt_layer_show_polygon.h"
#include "Qt_layer_show_greene_approximation.h"
#include "Qt_layer_show_ymonotone.h"
#include "Qt_layer_show_optimal_convex_partition.h"
#include "Qt_layer_show_polygon_points.h"

// icons
#include <CGAL/IO/pixmaps/ymonotone.xpm>
#include <CGAL/IO/pixmaps/greene_approx.xpm>
#include <CGAL/IO/pixmaps/show_polygon.xpm>
#include <CGAL/IO/pixmaps/optimal_convex.xpm>
#include <CGAL/IO/pixmaps/points.xpm>


namespace CGAL {
  Layers_toolbar::Layers_toolbar(Qt_widget *w, QMainWindow *mw, Cgal_Polygon *p) : 
     nr_of_buttons(0)
  {
    showP = new Qt_layer_show_polygon<Cgal_Polygon>(*p);
    showGA = new Qt_layer_show_greene_approx<Cgal_Polygon>(*p);
    showYM = new Qt_layer_show_ymonotone<Cgal_Polygon>(*p);
    showOC = new Qt_layer_show_optimal_convex<Cgal_Polygon>(*p);
    showPP = new Qt_layer_show_polygon_points<Cgal_Polygon>(*p);

    //set the widget
    widget = w;
    window = mw;

    widget->attach(showP);
    widget->attach(showGA);
    widget->attach(showYM);
    widget->attach(showOC);
    widget->attach(showPP);

    showGA->deactivate();
    showYM->deactivate();
    showOC->deactivate();

    maintoolbar = new QToolBar("tools", mw, QMainWindow::Top, TRUE, "Tools");
		

    but[0] = new QToolButton(maintoolbar, "show_polygon");
    but[0]->setPixmap(QPixmap( (const char**)show_polygon_xpm ));
    but[0]->setTextLabel("Show Simple Polygon");
    but[1] = new QToolButton(maintoolbar, "greene_approx");
    but[1]->setTextLabel("Show Greene Approximation");
    but[1]->setPixmap(QPixmap( (const char**)greene_approx_xpm ));
    but[2] = new QToolButton(maintoolbar, "ymonotone");
    but[2]->setTextLabel("Show Y Monotone Partition");
    but[2]->setPixmap(QPixmap( (const char**)ymonotone_xpm ));
    but[3] = new QToolButton(maintoolbar, "optimal_convex");
    but[3]->setPixmap(QPixmap( (const char**)optimal_convex_xpm ));
    but[3]->setTextLabel("Show Optimal Convex Partition");
    but[4] = new QToolButton(maintoolbar, "show_points");
    but[4]->setPixmap(QPixmap( (const char**)points_xpm ));
    but[4]->setTextLabel("Show Polygon Vertices");
    
    nr_of_buttons = 5;
    button_group = new QButtonGroup(0, "nonexclusive");
    
    for(int i =0; i<nr_of_buttons; i++){
      but[i]->setToggleButton(true);
      but[i]->toggle();
      button_group->insert(but[i]);
    }
    but[1]->toggle();
    but[2]->toggle();
    but[3]->toggle();
    connect(but[0], SIGNAL(stateChanged(int)),
        showP, SLOT(stateChanged(int)));
    connect(but[1], SIGNAL(stateChanged(int)),
        showGA, SLOT(stateChanged(int)));
    connect(but[2], SIGNAL(stateChanged(int)),
        showYM, SLOT(stateChanged(int)));
    connect(but[3], SIGNAL(stateChanged(int)),
        showOC, SLOT(stateChanged(int)));
    connect(but[4], SIGNAL(stateChanged(int)),
        showPP, SLOT(stateChanged(int)));
    connect(button_group, SIGNAL(clicked(int)),
          widget, SLOT(redraw()));
  }	

  Layers_toolbar::~Layers_toolbar()
  {
    delete showP;
    delete showGA;
    delete showYM;
    delete showOC;
    delete showPP;
    delete button_group;
  };

}//end namespace

#include "Qt_widget_toolbar_layers.moc"

#endif
