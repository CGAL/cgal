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
// file          : src/Qt_widget_toolbar_layers.C
// package       : Qt_widget
// author(s)     : Ursu Radu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


#ifdef CGAL_USE_QT

#include "Qt_widget_toolbar_layers.h"
#include <CGAL/IO/Qt_layer_show_mouse_coordinates.h>

#include "Qt_layer_show_polygon.h"
#include "Qt_layer_show_greene_approximation.h"
#include "Qt_layer_show_ymonotone.h"
#include "Qt_layer_show_optimal_convex_partition.h"
#include "Qt_layer_show_polygon_points.h"

// icons
#include <CGAL/IO/pixmaps/mouse_coord.xpm>
#include <CGAL/IO/pixmaps/ymonotone.xpm>
#include <CGAL/IO/pixmaps/greene_approx.xpm>
#include <CGAL/IO/pixmaps/show_polygon.xpm>
#include <CGAL/IO/pixmaps/optimal_convex.xpm>
#include <CGAL/IO/pixmaps/points.xpm>


namespace CGAL {
  Layers_toolbar::Layers_toolbar(Qt_widget *w, QMainWindow *mw, Cgal_Polygon *p) : 
     nr_of_buttons(0)
  {
    showMC = new Qt_layer_mouse_coordinates(*mw);
    showP = new Qt_layer_show_polygon<Cgal_Polygon>(*p);
    showGA = new Qt_layer_show_greene_approx<Cgal_Polygon>(*p);
    showYM = new Qt_layer_show_ymonotone<Cgal_Polygon>(*p);
    showOC = new Qt_layer_show_optimal_convex<Cgal_Polygon>(*p);
    showPP = new Qt_layer_show_polygon_points<Cgal_Polygon>(*p);

    //set the widget
    widget = w;
    window = mw;
    window->statusBar();

    widget->attach(showMC);
    widget->attach(showP);
    widget->attach(showGA);
    widget->attach(showYM);
    widget->attach(showOC);
    widget->attach(showPP);

    showGA->deactivate();
    showYM->deactivate();
    showOC->deactivate();

    maintoolbar = new QToolBar("tools", mw, QMainWindow::Top, TRUE, "Tools");
		

    but[0] = new QToolButton(maintoolbar, "mouse coordinates");
    but[0]->setPixmap(QPixmap( (const char**)mouse_coord_xpm ));
    but[0]->setTextLabel("Show Mouse Coordinates");
    but[1] = new QToolButton(maintoolbar, "show_polygon");
    but[1]->setPixmap(QPixmap( (const char**)show_polygon_xpm ));
    but[1]->setTextLabel("Show Simple Polygon");
    but[2] = new QToolButton(maintoolbar, "greene_approx");
    but[2]->setTextLabel("Show Greene Approximation");
    but[2]->setPixmap(QPixmap( (const char**)greene_approx_xpm ));
    but[3] = new QToolButton(maintoolbar, "ymonotone");
    but[3]->setTextLabel("Show Y Monotone Partition");
    but[3]->setPixmap(QPixmap( (const char**)ymonotone_xpm ));
    but[4] = new QToolButton(maintoolbar, "optimal_convex");
    but[4]->setPixmap(QPixmap( (const char**)optimal_convex_xpm ));
    but[4]->setTextLabel("Show Optimal Convex Partition");
    but[5] = new QToolButton(maintoolbar, "show_points");
    but[5]->setPixmap(QPixmap( (const char**)points_xpm ));
    but[5]->setTextLabel("Show Polygon Vertices");
    
    nr_of_buttons = 6;		
    button_group = new QButtonGroup(0, "nonexclusive");
    
    for(int i =0; i<nr_of_buttons; i++){
      but[i]->setToggleButton(true);
      but[i]->toggle();
      button_group->insert(but[i]);
    }
    but[2]->toggle();
    but[3]->toggle();
    but[4]->toggle();
    connect(but[0], SIGNAL(stateChanged(int)),
        showMC, SLOT(stateChanged(int)));
    connect(but[1], SIGNAL(stateChanged(int)),
        showP, SLOT(stateChanged(int)));
    connect(but[2], SIGNAL(stateChanged(int)),
        showGA, SLOT(stateChanged(int)));
    connect(but[3], SIGNAL(stateChanged(int)),
        showYM, SLOT(stateChanged(int)));
    connect(but[4], SIGNAL(stateChanged(int)),
        showOC, SLOT(stateChanged(int)));
    connect(but[5], SIGNAL(stateChanged(int)),
        showPP, SLOT(stateChanged(int)));
    connect(button_group, SIGNAL(clicked(int)),
          widget, SLOT(redraw()));
  }	

  Layers_toolbar::~Layers_toolbar()
  {
    delete showMC;
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
