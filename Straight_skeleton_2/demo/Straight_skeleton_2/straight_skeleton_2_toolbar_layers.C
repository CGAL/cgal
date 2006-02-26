// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//

#ifdef CGAL_USE_QT

#include "straight_skeleton_2_toolbar_layers.h"
#include "straight_skeleton_2_layers.h"

// icons
#include <CGAL/IO/pixmaps/voronoi.xpm>
#include <CGAL/IO/pixmaps/show_polygon.xpm>
#include <CGAL/IO/pixmaps/polygon.xpm>

#include <qiconset.h>

Layers_toolbar::Layers_toolbar(CGAL::Qt_widget*     w
                              ,QMainWindow*         mw
                              ,demo::Regions const& in
                              ,demo::Sls     const& sls
                              ,demo::Regions const& out
                              ) : QToolBar(mw, "LT"),
     nr_of_buttons(0)
  {
    showI   = new Qt_layer_show_regions <demo::Regions>(in,CGAL::RED);
    showSLS = new Qt_layer_show_skeleton<demo::Sls>    (sls);
    showO   = new Qt_layer_show_regions <demo::Regions>(out,CGAL::BLACK);

    //set the widget
    widget = w;
    window = mw;

    widget->attach(showI);
    widget->attach(showSLS);
    widget->attach(showO);


    QIconSet set0(QPixmap( (const char**)show_polygon_small_xpm ),
                  QPixmap( (const char**)show_polygon_xpm ));
    QIconSet set1(QPixmap( (const char**)voronoi_small_xpm ),
                  QPixmap( (const char**)voronoi_xpm ));
    QIconSet set2(QPixmap( (const char**)polygon_small_xpm ),
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

    nr_of_buttons = 3;
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
        showSLS, SLOT(stateChanged(int)));
    connect(but[2], SIGNAL(stateChanged(int)),
        showO, SLOT(stateChanged(int)));
  }

  Layers_toolbar::~Layers_toolbar()
  {
    delete showI;
    delete showSLS;
    delete showO;
    delete button_group;
  };


#include "straight_skeleton_2_toolbar_layers.moc"

#endif
