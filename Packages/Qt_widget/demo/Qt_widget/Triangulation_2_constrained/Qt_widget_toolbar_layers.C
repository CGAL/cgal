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
#include <CGAL/IO/pixmaps/constrained.xpm>
#include <CGAL/IO/pixmaps/triangulation.xpm>
#include <CGAL/IO/pixmaps/mouse_coord.xpm>


namespace CGAL {
  Layers_toolbar::Layers_toolbar(Qt_widget *w, QMainWindow *mw, CDT1 *t) : 
     nr_of_buttons(0)
  {
    showT   = new Qt_layer_show_triangulation< CDT1 >(*t);
    showP   = new Qt_layer_show_points< CDT1 >(*t);
    showMC  = new Qt_layer_mouse_coordinates(*mw);
    showC   = new Qt_layer_show_constraineds<CDT1>(*t);

    //set the widget
    widget = w;
    window = mw;
    window->statusBar();

    widget->attach(showT);    
    widget->attach(showC);
    widget->attach(showP);
    widget->attach(showMC);
    

    maintoolbar = new QToolBar("tools", mw, QMainWindow::Top, TRUE, "Tools");
		
    but[0] = new QToolButton(maintoolbar, "triangulation");
    but[0]->setPixmap(QPixmap( (const char**)triangulation_xpm ));
    but[0]->setTextLabel("Show Triangulation");
    but[1] = new QToolButton(maintoolbar, "constraineds");
    but[1]->setPixmap(QPixmap( (const char**)constrained_xpm ));
    but[1]->setTextLabel("Show Constraineds");
    but[2] = new QToolButton(maintoolbar, "vertices");
    but[2]->setPixmap(QPixmap( (const char**)points_xpm ));
    but[2]->setTextLabel("Show Vertices");
    but[3] = new QToolButton(maintoolbar, "mouse_coord");
    but[3]->setPixmap(QPixmap( (const char**)mouse_coord_xpm ));
    but[3]->setTextLabel("Mouse Coordinates");
		

    nr_of_buttons = 4;
	  button_group = new QButtonGroup(0, "nonexclusive");
    for(int i =0; i<nr_of_buttons; i++)
    {
      but[i]->setToggleButton(TRUE);
      button_group->insert(but[i]);
      but[i]->toggle();
    }    
    connect(button_group, SIGNAL(clicked(int)),
          widget, SLOT(redraw()));
    
    connect(but[0], SIGNAL(stateChanged(int)),
        showT, SLOT(stateChanged(int)));
    connect(but[1], SIGNAL(stateChanged(int)),
        showC, SLOT(stateChanged(int)));
    connect(but[2], SIGNAL(stateChanged(int)),
        showP, SLOT(stateChanged(int)));
    connect(but[3], SIGNAL(stateChanged(int)),
        showMC, SLOT(stateChanged(int)));
  }
  
}//end namespace

#include "Qt_widget_toolbar_layers.moc"

#endif
