// Copyright (c) 2002  ETH Zurich (Switzerland).
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
#include <CGAL/IO/pixmaps/min_rectangle.xpm>
#include <CGAL/IO/pixmaps/min_parallelogram.xpm>

#include <qiconset.h>

/* XPM */
static const char * lines_small_xpm[] = {
"16 16 3 1",
" 	c None",
".	c #35E8D9",
"+	c #000000",
"    ........    ",
"   ..+.......   ",
"  ....+.......  ",
" ......+....... ",
"........+.......",
".+.......+......",
"..+.......+.....",
"...+.......+....",
"....+.......+...",
".....+.......+..",
"......+.......+.",
".......+........",
" .......+...... ",
"  .......+....  ",
"   .......+..   ",
"    ........    "};



Layers_toolbar::Layers_toolbar(CGAL::Qt_widget *w, QMainWindow *mw,
			       std::list<Point_2> *l_of_p) :
  QToolBar(mw, "LT"), widget(w), window(mw), nr_of_buttons(0)
  {

    showPL  = new Qt_layer_show_parallelogram<Rp>(l_of_p);
    showP   = new Qt_layer_show_points<Rp>(l_of_p);
    showLS  = new Qt_layer_show_strip<Rp>(l_of_p);
    showR   = new Qt_layer_show_rectangle<Rp>(l_of_p);

    //set the widget

    widget->attach(showR);
    widget->attach(showPL);
    widget->attach(showP);
    widget->attach(showLS);

    QIconSet set0(QPixmap( (const char**)points_small_xpm ),
                  QPixmap( (const char**)points_xpm ));
    QIconSet set1(QPixmap( (const char**)min_parallelogram_small_xpm ),
                  QPixmap( (const char**)min_parallelogram_xpm ));
    QIconSet set2(QPixmap( (const char**)lines_small_xpm ),
                  QPixmap( (const char**)lines_small_xpm ));
    QIconSet set3(QPixmap( (const char**)min_rectangle_small_xpm ),
                  QPixmap( (const char**)min_rectangle_xpm ));

    but[0] = new QToolButton(this, "points");
    but[0]->setIconSet(set0);
    but[0]->setTextLabel("Show Points");
    but[1] = new QToolButton(this, "Minimum_parallelogram");
    but[1]->setIconSet(set1);
    but[1]->setTextLabel("Show Minimum Parallelogram");
    but[2] = new QToolButton(this, "Show_lines");
    but[2]->setIconSet(set2);
    but[2]->setTextLabel("Show Lines");
    but[3] = new QToolButton(this, "Minimum_rectangle");
    but[3]->setIconSet(set3);
    but[3]->setTextLabel("Show Minimum Rectangle");

    nr_of_buttons = 4;
    button_group = new QButtonGroup(0, "nonexclusive");
    for(int i =0; i<nr_of_buttons; i++)
    {
	but[i]->setToggleButton(TRUE);
	but[i]->toggle();
	button_group->insert(but[i]);
    }
    connect(but[0], SIGNAL(stateChanged(int)),
      showP, SLOT(stateChanged(int)));
    connect(but[1], SIGNAL(stateChanged(int)),
      showPL, SLOT(stateChanged(int)));
    connect(but[2], SIGNAL(stateChanged(int)),
      showLS, SLOT(stateChanged(int)));
    connect(but[3], SIGNAL(stateChanged(int)),
      showR, SLOT(stateChanged(int)));

    connect(button_group, SIGNAL(clicked(int)),
      widget, SLOT(redraw()));

  }//end constructor

#include "Qt_widget_toolbar_layers.moc"

