// Copyright (c) 2002  INRIA Sophia-Antipolis (France).
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

#include "Qt_widget_toolbar_layers.h"

// icons
#include <CGAL/IO/pixmaps/points.xpm>
#include <CGAL/IO/pixmaps/voronoi.xpm>
#include <CGAL/IO/pixmaps/triangulation.xpm>
#include <CGAL/IO/pixmaps/alpha_shape.xpm>

#include <qiconset.h>


/* XPM */
static const char * image_small_xpm[] = {
"16 16 6 1",
"       c None",
".      c #35E8D9",
"+      c #000000",
"@      c #E22F2F",
"#      c #F7F309",
"$      c #0C4DF2",
"    ..++....    ",
"   ...++.....   ",
"  ............  ",
" ..@@@++@...... ",
"..@@@@++@@@.....",
"..@@@@++@@@@....",
"...@@@++@@@@@...",
"..####++###@@...",
"...###++#####...",
"..$$$#++#####...",
"..$$$$++######..",
"...$$$++$$$$$...",
" ...$$$+++$$$.. ",
"  ...$$$++$$..  ",
"   ...$$$$$..   ",
"    ........    "};



/* XPM */
static const char * image_xpm[] = {
"32 32 5 1",
"       c None",
".      c #000000",
"+      c #F23818",
"@      c #F9ED09",
"#      c #0C69F4",
"                                ",
"                                ",
"           .......              ",
"           .......              ",
"           .......      +++++   ",
"         +++.....++   ++++++++  ",
"    +++++++++++++++++++++++++++ ",
"   ++++++++++++++++++++++++++++ ",
"  +++++++++++....++++++++++++++ ",
"  +++++++++++....++++++++++++++ ",
"  +++++++++++....++++++++++++++ ",
"  +++++++++++....+++++++++++++  ",
"  ++++++@@@@@....++++++++++++   ",
"   @@@@@@@@@@....@++++@@@@@+    ",
"  @@@@@@@@@@@....@@++@@@@@@@    ",
" @@@@@@@@@@@@....@@@@@@@@@@@@   ",
" @@@@@@@@@@@@....@@@@@@@@@@@@   ",
"@@@@@@@@@@@@@....@@@@@@@@@@@@   ",
"@@@@@@@@@@@@@....@@@@@@@@@@@@   ",
"@@@@@#######@....@@@@@@@@@@@@   ",
"@@@@#########....@@@@@@@@@@@    ",
"@@@##########....#@@#########   ",
" @@##########....#############  ",
"  ###########....############## ",
"  ###########.....############# ",
"  ############.....############ ",
"  ############......########### ",
"  #############.....########### ",
"   ##############...##########  ",
"    ##### ###################   ",
"               #########        ",
"                                "};

Layers_toolbar::Layers_toolbar(CGAL::Qt_widget *w, QMainWindow *mw,
                                Delaunay *t, Alpha_shape *a, QImage
			       *i) : QToolBar(mw, "LT") ,nr_of_buttons(0)
{

  showT = new Qt_layer_show_triangulation< Delaunay >(*t);
  showV = new Qt_layer_show_voronoi< Delaunay >(*t);
  showP = new Qt_layer_show_points< Delaunay >(*t);
  showA = new Qt_layer_show_alpha_shape(a);
  showI = new Qt_layer_show_image(i);

  //set the widget
  widget = w;
  window = mw;

  widget->attach(showI);
  widget->attach(showT);
  widget->attach(showA);
  widget->attach(showV);
  widget->attach(showP);

  showV->deactivate();

  QIconSet set0(QPixmap( (const char**)image_small_xpm ),
                QPixmap( (const char**)image_xpm ));
  QIconSet set1(QPixmap( (const char**)triangulation_small_xpm ),
                QPixmap( (const char**)triangulation_xpm ));
  QIconSet set2(QPixmap( (const char**)alpha_shape_small_xpm ),
                QPixmap( (const char**)alpha_shape_xpm ));
  QIconSet set3(QPixmap( (const char**)voronoi_small_xpm ),
                QPixmap( (const char**)voronoi_xpm ));
  QIconSet set4(QPixmap( (const char**)points_small_xpm ),
                QPixmap( (const char**)points_xpm ));

  but[0] = new QToolButton(this, "image layer");
  but[0]->setIconSet(set0);
  but[0]->setTextLabel("Show Image");
  but[1] = new QToolButton(this, "triangulation layer");
  but[1]->setIconSet(set1);
  but[1]->setTextLabel("Show Triangulation");
  but[2] = new QToolButton(this, "alpha_shape layer");
  but[2]->setIconSet(set2);
  but[2]->setTextLabel("Show Alpha Shape");
  but[3] = new QToolButton(this, "voronoi layer");
  but[3]->setIconSet(set3);
  but[3]->setTextLabel("Show Voronoi");
  but[4] = new QToolButton(this, "vertices layer");
  but[4]->setIconSet(set4);
  but[4]->setTextLabel("Show Vertices");

  nr_of_buttons = 5;
  button_group = new QButtonGroup(0, "nonexclusive");
  for(int i =0; i<nr_of_buttons; i++)
  {
    but[i]->setToggleButton(TRUE);
    but[i]->toggle();
    button_group->insert(but[i]);
  }
  but[3]->toggle();
  connect(button_group, SIGNAL(clicked(int)),
        widget, SLOT(redraw()));

  connect(but[0], SIGNAL(stateChanged(int)),
      showI, SLOT(stateChanged(int)));
  connect(but[1], SIGNAL(stateChanged(int)),
      showT, SLOT(stateChanged(int)));
  connect(but[2], SIGNAL(stateChanged(int)),
      showA, SLOT(stateChanged(int)));
  connect(but[3], SIGNAL(stateChanged(int)),
      showV, SLOT(stateChanged(int)));
  connect(but[4], SIGNAL(stateChanged(int)),
      showP, SLOT(stateChanged(int)));

}

#include "Qt_widget_toolbar_layers.moc"

