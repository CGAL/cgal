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
// file          : triangulation_2_toolbar_layers.C
// package       : Qt_widget
// author(s)     : Ursu Radu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


#ifdef CGAL_USE_QT

#include "regular_triangulation_2_toolbar_layers.h"

// icons
#include <CGAL/IO/pixmaps/points.xpm>
#include <CGAL/IO/pixmaps/voronoi.xpm>
#include <CGAL/IO/pixmaps/triangulation.xpm>

#include <qiconset.h>


  Layers_toolbar::Layers_toolbar(CGAL::Qt_widget *w, QMainWindow *mw, 
                                 Regular_triangulation *t) :
    QToolBar(mw, "LT"), dt(t), nr_of_buttons(0)
  {
    showT   = new Qt_layer_show_triangulation< Regular_triangulation >(*t);
    showV   = new Qt_layer_show_voronoi< Regular_triangulation >(*t);
    showP   = new Qt_layer_show_points< Regular_triangulation >(*t);

    //set the widget
    widget = w;
    window = mw;

    widget->attach(showT);
    widget->attach(showV);
    widget->attach(showP);

    QIconSet set0(QPixmap( (const char**)triangulation_small_xpm ),
                  QPixmap( (const char**)triangulation_xpm ));
    QIconSet set1(QPixmap( (const char**)voronoi_small_xpm ),
                  QPixmap( (const char**)voronoi_xpm ));
    QIconSet set3(QPixmap( (const char**)points_small_xpm ),
                  QPixmap( (const char**)points_xpm ));

    but[0] = new QToolButton(this, "triangulation");
    but[0]->setIconSet(set0);
    but[0]->setTextLabel("Triangulation");
    but[1] = new QToolButton(this, "voronoi");
    but[1]->setIconSet(set1);
    but[1]->setTextLabel("Voronoi Diagram");
    but[2] = new QToolButton(this, "vertices");
    but[2]->setIconSet(set3);
    but[2]->setTextLabel("Vertices");

    nr_of_buttons = 3;
	  button_group = new QButtonGroup(0, "nonexclusive");
    for(int i =0; i<nr_of_buttons; i++)
    {
      but[i]->setToggleButton(true);
      but[i]->toggle();
      button_group->insert(but[i]);
    }
    connect(button_group, SIGNAL(clicked(int)),
          widget, SLOT(redraw()));
    
    connect(but[0], SIGNAL(stateChanged(int)),
        showT, SLOT(stateChanged(int)));
    connect(but[1], SIGNAL(stateChanged(int)),
        showV, SLOT(stateChanged(int)));
    connect(but[2], SIGNAL(stateChanged(int)),
        showP, SLOT(stateChanged(int)));
  }
  

#include "regular_triangulation_2_toolbar_layers.moc"

#endif
