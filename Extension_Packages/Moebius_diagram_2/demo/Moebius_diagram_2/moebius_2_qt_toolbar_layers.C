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

#include "moebius_2_qt_toolbar_layers.h"

// icons
#include <CGAL/IO/pixmaps/points.xpm>
#include <CGAL/IO/pixmaps/voronoi.xpm>
#include <CGAL/IO/pixmaps/triangulation.xpm>

#include <qiconset.h>


  Layers_toolbar::Layers_toolbar(CGAL::Qt_widget *w, QMainWindow *mw, 
                                 MD *d) :
    QToolBar(mw, "LT"), md(*d), nr_of_buttons(0)
  {
    //    showT   = new Qt_layer_show_triangulation< Regular_triangulation >(*t);
    showV   = new Qt_layer_show_moebius< MD >(d);
    showP   = new Qt_layer_show_points< MD> (d);
    showL   = new Qt_layer_show_lambdas< MD> (d);
    showM   = new Qt_layer_show_mus< MD> (d);

    //set the widget
    widget = w;
    window = mw;

    //    widget->attach(showT);
    widget->attach(showV);
    widget->attach(showP);
    widget->attach(showL);
    widget->attach(showM);

    QIconSet set0(QPixmap( (const char**)triangulation_small_xpm ),
                  QPixmap( (const char**)triangulation_xpm ));
    QIconSet set1(QPixmap( (const char**)voronoi_small_xpm ),
                  QPixmap( (const char**)voronoi_xpm ));
    QIconSet set3(QPixmap( (const char**)points_small_xpm ),
                  QPixmap( (const char**)points_xpm ));

    but[0] = new QToolButton(this, "moebius");
    but[0]->setIconSet(set1);
    but[0]->setTextLabel("Möbius Diagram");

    but[1] = new QToolButton(this, "vertices");
    but[1]->setIconSet(set3);
    but[1]->setTextLabel("Vertices");

    but[2] = new QToolButton(this, "lambda");
    but[2]->setIconSet(set3);
    but[2]->setTextLabel("1 / lambda");

    but[3] = new QToolButton(this, "mu");
    but[3]->setIconSet(set3);
    but[3]->setTextLabel("mu");

    nr_of_buttons = 4;
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
        showV, SLOT(stateChanged(int)));
    connect(but[1], SIGNAL(stateChanged(int)),
        showP, SLOT(stateChanged(int)));
    connect(but[2], SIGNAL(stateChanged(int)),
        showL, SLOT(stateChanged(int)));
    connect(but[3], SIGNAL(stateChanged(int)),
        showM, SLOT(stateChanged(int)));
  }
  

#include "moebius_2_qt_toolbar_layers.moc"

#endif
