// Copyright (c) 1997-2002  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Radu Ursu


#include <CGAL/basic.h>


#include "delaunay_triangulation_2_toolbar_layers.h"

// icons
#include <CGAL/IO/pixmaps/points.xpm>
#include <CGAL/IO/pixmaps/nearest_vertex.xpm>
#include <CGAL/IO/pixmaps/voronoi.xpm>
#include <CGAL/IO/pixmaps/triangulation.xpm>

#include <qiconset.h>

/* XPM */
static const char *circum_circle_xpm[] = {
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

/* XPM */
static const char * circum_circle_small_xpm[] = {
"16 16 4 1",
" 	c None",
".	c #35E8D9",
"+	c #F4072E",
"@	c #000000",
"    ........    ",
"   ..++++....   ",
"  ..+...@+@...  ",
" ..+.@@@..+@@.. ",
"..+@@.....@+.@..",
".@+......@.+..@.",
".@+@@....@.+..@.",
"..+..@..@..+@@@.",
"..@+..@.@.+...@.",
"...@+..@@+...@..",
"....@++++...@...",
"....@..@...@....",
" ...@..@..@.... ",
"  ...@.@.@....  ",
"   ..@.@@....   ",
"    ..@@....    "};

  Layers_toolbar::Layers_toolbar(CGAL::Qt_widget *w, QMainWindow *mw,
                                 Delaunay *t) : QToolBar(mw, "LT"),
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

    QIconSet set0(QPixmap( (const char**)triangulation_small_xpm ),
                  QPixmap( (const char**)triangulation_xpm ));
    QIconSet set1(QPixmap( (const char**)voronoi_small_xpm ),
                  QPixmap( (const char**)voronoi_xpm ));
    QIconSet set2(QPixmap( (const char**)nearest_vertex_small_xpm ),
                  QPixmap( (const char**)nearest_vertex_xpm ));
    QIconSet set3(QPixmap( (const char**)points_small_xpm ),
                  QPixmap( (const char**)points_xpm ));
    QIconSet set4(QPixmap( (const char**)circum_circle_small_xpm ),
                  QPixmap( (const char**)circum_circle_xpm ));

    but[0] = new QToolButton(this, "triangulation");
    but[0]->setIconSet(set0);
    but[0]->setTextLabel("Triangulation");
    but[1] = new QToolButton(this, "voronoi");
    but[1]->setIconSet(set1);
    but[1]->setTextLabel("Voronoi Diagram");
    but[2] = new QToolButton(this, "nearest_vertex");
    but[2]->setIconSet(set2);
    but[2]->setTextLabel("Nearest Vertex");
    but[3] = new QToolButton(this, "vertices");
    but[3]->setIconSet(set3);
    but[3]->setTextLabel("Vertices");
    but[4] = new QToolButton(this, "circles");
    but[4]->setIconSet(set4);
    but[4]->setTextLabel("Circumscribed Circle");

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


#include "delaunay_triangulation_2_toolbar_layers.moc"

