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
#include <CGAL/IO/pixmaps/line.xpm>
#include <CGAL/IO/pixmaps/min_rectangle.xpm>
#include <CGAL/IO/pixmaps/min_parallelogram.xpm>
#include <CGAL/IO/pixmaps/mouse_coord.xpm>


namespace CGAL {
  Layers_toolbar::Layers_toolbar(Qt_widget *w, QMainWindow *mw, std::list<Point>	*l_of_p) : 
    widget(w), window(mw), nr_of_buttons(0)
  {
      
    showMC  = new Qt_layer_mouse_coordinates(*mw);
    showPL  = new Qt_layer_show_parallelogram<Rp>(l_of_p);
    showP   = new Qt_layer_show_points<Rp>(l_of_p);
    showLS  = new Qt_layer_show_strip<Rp>(l_of_p);
    showR   = new Qt_layer_show_rectangle<Rp>(l_of_p);
    
    //set the widget
    window->statusBar();

    widget->attach(showMC);
    widget->attach(showR);
    widget->attach(showPL);
    widget->attach(showP);
    widget->attach(showLS);

    maintoolbar = new QToolBar("tools", mw, QMainWindow::Top, TRUE, "Tools");

    but[0] = new QToolButton(maintoolbar, "mouse_coord");
    but[0]->setPixmap(QPixmap( (const char**)mouse_coord_xpm ));
    but[1] = new QToolButton(maintoolbar, "points");
    but[1]->setPixmap(QPixmap( (const char**)points_xpm ));
    but[2] = new QToolButton(maintoolbar, "Minimum_parallelogram");
    but[2]->setPixmap(QPixmap( (const char**)min_parallelogram_xpm ));
    but[3] = new QToolButton(maintoolbar, "Show_lines");
    but[3]->setPixmap(QPixmap( (const char**)line_xpm ));
    but[4] = new QToolButton(maintoolbar, "Minimum_rectangle");
    but[4]->setPixmap(QPixmap( (const char**)min_rectangle_xpm ));

    nr_of_buttons = 5;
    button_group = new QButtonGroup(0, "nonexclusive");	
    for(int i =0; i<nr_of_buttons; i++)
    {
	but[i]->setToggleButton(TRUE);
	but[i]->toggle();
	button_group->insert(but[i]);
    }
    connect(but[0], SIGNAL(stateChanged(int)), 
      showMC, SLOT(stateChanged(int)));
    connect(but[1], SIGNAL(stateChanged(int)),
      showP, SLOT(stateChanged(int)));
    connect(but[2], SIGNAL(stateChanged(int)), 
      showPL, SLOT(stateChanged(int)));
    connect(but[3], SIGNAL(stateChanged(int)), 
      showLS, SLOT(stateChanged(int)));
    connect(but[4], SIGNAL(stateChanged(int)), 
      showR, SLOT(stateChanged(int)));

    connect(button_group, SIGNAL(clicked(int)),
      widget, SLOT(redraw()));
    
  }//end constructor

}//end namespace

#include "Qt_widget_toolbar_layers.moc"

#endif
