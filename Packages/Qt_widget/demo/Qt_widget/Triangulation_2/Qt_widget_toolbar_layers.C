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

/* XPM */
static char *circum_circle_xpm[] = {
/* columns rows colors chars-per-pixel */
"32 32 5 1",
"  c opaque",
". c green",
"X c red",
"o c #c0c0c0",
"O c None",
/* pixels */
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOXXXXXOOOOOOOOOOOOOOO",
"OOOOOOOOOXXXOOOOOXXXOOOOOOOOOOOO",
"OOOOOOOXXOOOOOOOOOOOXXOOOOOOOOOO",
"OOOOOOXOOOOOOOOOOOOOOOXOOOOOOOOO",
"OOOOOXOOOOOOOOOOOOOOOOOXOOOOOOOO",
"OOOOXOOOOOOOOOOOOOOOOOOOXOOOOOOO",
"OOOXOOOOOOOOOOOOOOOOOOOOOXOOOOOO",
"OOOXOOOOOOOOOOOOOOOOOOOOOXOOOOOO",
"OOOXOOOOOOOOOOOOOOOOOOOOOXOOOOOO",
"OOXOOOOOOOOOOOOOOOOOOOOOOOXOOOOO",
"OOXOOOOOOOOOOOOOOOOOOOOOOOXOOOOO",
"OOXOOOOOOOOOOOOOOOOOOOOOOOXOOOOO",
"OO..                     ..OOOOO",
"OO.. ooooooooooooooooooo ..OOOOO",
"OO Xo ooooooooooooooooo oX OOOOO",
"OO Xoo  ooooooooooooo  ooX OOOOO",
"OO oXooo ooooooooooo oooXo OOOOO",
"OO ooXooo  ooooooo  oooXoo OOOOO",
"OO oooXoooo ooooo ooooXooo OOOOO",
"OO ooooXXooo  o  oooXXoooo OOOOO",
"OO ooooooXXXo..ooXXXoooooo OOOOO",
"OO oooooooooX..XXooooooooo OOOOO",
"OO ooooooo  ooooo oooooooo OOOOO",
"OO oooooo oooooooo  oooooo OOOOO",
"OO oooo  ooooooooooo ooooo OOOOO",
"OO ooo oooooooooooooo  ooo OOOOO",
"OO o  ooooooooooooooooo oo OOOOO",
"OO  oooooooooooooooooooo   OOOOO",
"OO                         OOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"
};

namespace CGAL {
  Layers_toolbar::Layers_toolbar(Qt_widget *w, QMainWindow *mw, Delaunay *t) : 
    dt(t), nr_of_buttons(0)
  {
    showT   = new Qt_layer_show_triangulation< Delaunay >(*t);
    showV   = new Qt_layer_show_voronoi< Delaunay >(*t);
    showP   = new Qt_layer_show_points< Delaunay >(*t);
    showNV  = new Qt_layer_nearest_vertex< Delaunay >(*t);
    showCC  = new Qt_layer_circum_circle< Delaunay >(*t);

    //set the widget
    widget = w;
    window = mw;

    widget->attach(showT);
    widget->attach(showV);
    widget->attach(showNV);
    widget->attach(showP);
    widget->attach(showCC);
    showNV->deactivate();
    showCC->deactivate();

    maintoolbar = new QToolBar("tools", mw, QMainWindow::Top, TRUE, "Tools");
		
    but[0] = new QToolButton(maintoolbar, "triangulation");
    but[0]->setPixmap(QPixmap( (const char**)triangulation_xpm ));
    but[0]->setTextLabel("Triangulation");
    but[1] = new QToolButton(maintoolbar, "voronoi");
    but[1]->setPixmap(QPixmap( (const char**)voronoi_xpm ));
    but[1]->setTextLabel("Voronoi Diagram");
    but[2] = new QToolButton(maintoolbar, "nearest_vertex");
    but[2]->setPixmap(QPixmap( (const char**)nearest_vertex_xpm ));
    but[2]->setTextLabel("Nearest Vertex");
    but[3] = new QToolButton(maintoolbar, "vertices");
    but[3]->setPixmap(QPixmap( (const char**)points_xpm ));
    but[3]->setTextLabel("Vertices");
    but[4] = new QToolButton(maintoolbar, "circles");
    but[4]->setPixmap(QPixmap( (const char**)circum_circle_xpm ));
    but[4]->setTextLabel("Circuscribed Circle");
		
    nr_of_buttons = 5;
	  button_group = new QButtonGroup(0, "nonexclusive");
    for(int i =0; i<nr_of_buttons; i++)
    {
      but[i]->setToggleButton(TRUE);
      but[i]->toggle();
      button_group->insert(but[i]);
    }
    but[2]->toggle();
    but[4]->toggle();
    connect(button_group, SIGNAL(clicked(int)),
          widget, SLOT(redraw()));
    
    connect(but[0], SIGNAL(stateChanged(int)),
        showT, SLOT(stateChanged(int)));
    connect(but[1], SIGNAL(stateChanged(int)),
        showV, SLOT(stateChanged(int)));
    connect(but[2], SIGNAL(stateChanged(int)),
        showNV, SLOT(stateChanged(int)));
    connect(but[3], SIGNAL(stateChanged(int)),
        showP, SLOT(stateChanged(int)));
    connect(but[4], SIGNAL(stateChanged(int)),
        showCC, SLOT(stateChanged(int)));

  }
  
}//end namespace

#include "Qt_widget_toolbar_layers.moc"

#endif
