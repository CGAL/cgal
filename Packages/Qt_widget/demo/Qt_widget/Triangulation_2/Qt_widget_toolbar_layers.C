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

// icons
#include <CGAL/IO/pixmaps/points.xpm>
#include <CGAL/IO/pixmaps/nearest_vertex.xpm>
#include <CGAL/IO/pixmaps/voronoi.xpm>
#include <CGAL/IO/pixmaps/triangulation.xpm>
#include <CGAL/IO/pixmaps/mouse_coord.xpm>


namespace CGAL {
  Layers_toolbar::Layers_toolbar(Qt_widget *w, QMainWindow *mw, Delaunay *t) : 
    dt(t), nr_of_buttons(0)
  {
    showT   = new Qt_layer_show_triangulation< Delaunay >(*t);
    showV   = new Qt_layer_show_voronoi< Delaunay >(*t);
    showP   = new Qt_layer_show_points< Delaunay >(*t);
    showNV  = new Qt_layer_nearest_vertex< Delaunay >(*t);
    showMC  = new Qt_layer_mouse_coordinates(*mw);

    //set the widget
    widget = w;
    window = mw;
    window->statusBar();

    widget->attach(showT);
    widget->attach(showV);
    widget->attach(showNV);
    widget->attach(showP);
    widget->attach(showMC);
    showNV->deactivate();
    

    maintoolbar = new QToolBar("tools", mw, QMainWindow::Top, TRUE, "Tools");
		
    but[0] = new QToolButton(maintoolbar, "triangulation");
    but[0]->setPixmap(QPixmap( (const char**)triangulation_xpm ));
    but[1] = new QToolButton(maintoolbar, "voronoi");
    but[1]->setPixmap(QPixmap( (const char**)voronoi_xpm ));
    but[2] = new QToolButton(maintoolbar, "nearest_vertex");
    but[2]->setPixmap(QPixmap( (const char**)nearest_vertex_xpm ));
    but[3] = new QToolButton(maintoolbar, "vertices");
    but[3]->setPixmap(QPixmap( (const char**)points_xpm ));
    but[4] = new QToolButton(maintoolbar, "mouse_coord");
    but[4]->setPixmap(QPixmap( (const char**)mouse_coord_xpm ));
		
    nr_of_buttons = 5;
	  button_group = new QButtonGroup(0, "nonexclusive");
    for(int i =0; i<nr_of_buttons; i++)
    {
      but[i]->setToggleButton(TRUE);
      but[i]->toggle();
      button_group->insert(but[i]);
    }
    but[2]->toggle();
    connect(button_group, SIGNAL(clicked(int)),
          this, SLOT(redraw_win(int)));
    
    connect(but[0], SIGNAL(stateChanged(int)),
        showT, SLOT(stateChanged(int)));
    connect(but[1], SIGNAL(stateChanged(int)),
        showV, SLOT(stateChanged(int)));
    connect(but[2], SIGNAL(stateChanged(int)),
        showNV, SLOT(stateChanged(int)));
    connect(but[3], SIGNAL(stateChanged(int)),
        showP, SLOT(stateChanged(int)));
    connect(but[4], SIGNAL(stateChanged(int)),
        showMC, SLOT(stateChanged(int)));



  }
  void Layers_toolbar::redraw_win(int i){
    widget->redraw();
  }

  
}//end namespace

#include "Qt_widget_toolbar_layers.moc"

#endif
