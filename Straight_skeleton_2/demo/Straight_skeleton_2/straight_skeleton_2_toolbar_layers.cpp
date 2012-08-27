// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//

#include <CGAL/basic.h>


#include "cgal_types.h"
#include "ss_types.h"
#include "straight_skeleton_2_toolbar_layers.h"
#include "straight_skeleton_2_layers.h"

// icons
#include <CGAL/IO/pixmaps/voronoi.xpm>
#include <CGAL/IO/pixmaps/show_polygon.xpm>
#include <CGAL/IO/pixmaps/polygon.xpm>

#include <qiconset.h>

Layers_toolbar::Layers_toolbar(CGAL::Qt_widget*        w
                              ,QMainWindow*            mw
                              ,demo::Regions  const&   in
                              ,demo::SSkelPtr const&   sskel
                              ,demo::Regions  const&   out
                              ) : QToolBar(mw, "LT"),
     nr_of_buttons(0)
  {
    showI     = new Qt_layer_show_regions <demo::Regions> (this,"Input",in,CGAL::RED);
    showSSkel = new Qt_layer_show_skeleton<demo::SSkel>   (this,"Skeleton",sskel);
    showO     = new Qt_layer_show_regions <demo::Regions> (this,"Offset",out,CGAL::BLACK);
    progress  = new Qt_layer_show_progress                (this,"Progress");

    //set the widget
    widget = w;
    window = mw;

    widget->attach(showI);
    widget->attach(progress);
    widget->attach(showSSkel);
    widget->attach(showO);

    QIconSet set0(QPixmap( (const char**)show_polygon_small_xpm ),
                  QPixmap( (const char**)show_polygon_xpm ));
    QIconSet set1(QPixmap( (const char**)voronoi_small_xpm ),
                  QPixmap( (const char**)voronoi_xpm ));
    QIconSet set2(QPixmap( (const char**)polygon_small_xpm ),
                  QPixmap( (const char**)polygon_xpm ));
    QIconSet set3(QPixmap( (const char**)polygon_small_xpm ),
                  QPixmap( (const char**)polygon_xpm ));

    but[0] = new QToolButton(this, "polygon");
    but[0]->setIconSet(set0);
    but[0]->setTextLabel("Show Simple Polygon");

    but[1] = new QToolButton(this, "straight_skeleton");
    but[1]->setIconSet(set1);
    but[1]->setTextLabel("Show Straight Skeleton");

    but[2] = new QToolButton(this, "offset");
    but[2]->setIconSet(set2);
    but[2]->setTextLabel("Show Polygon Offset");

    but[3] = new QToolButton(this, "progress");
    but[3]->setIconSet(set3);
    but[3]->setTextLabel("Show Progress");

    nr_of_buttons = 4;
    button_group = new QButtonGroup(0, "nonexclusive");

    for(int i =0; i<nr_of_buttons; i++){
      but[i]->setToggleButton(true);
      but[i]->toggle();
      button_group->insert(but[i]);
    }
    //but[1]->toggle();
    connect(but[0], SIGNAL(stateChanged(int)),
        showI, SLOT(stateChanged(int)));
    connect(but[1], SIGNAL(stateChanged(int)),
        showSSkel, SLOT(stateChanged(int)));
    connect(but[2], SIGNAL(stateChanged(int)),
        showO, SLOT(stateChanged(int)));
    connect(but[3], SIGNAL(stateChanged(int)),
        progress, SLOT(stateChanged(int)));
  }

  Layers_toolbar::~Layers_toolbar()
  {
    delete showI;
    delete showSSkel;
    delete showO;
    delete progress;
    delete button_group;
  }


#include "straight_skeleton_2_toolbar_layers.moc"

