// Copyright (c) 2002  Max Planck Institut fuer Informatik (Germany).
// All rights reserved.
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
//
// Author(s)     : Radu Ursu

#include <CGAL/basic.h>


#include "partition_2_toolbar_layers.h"
#include "partition_2_layers.h"

// icons
#include <CGAL/IO/pixmaps/ymonotone.xpm>
#include <CGAL/IO/pixmaps/greene_approx.xpm>
#include <CGAL/IO/pixmaps/show_polygon.xpm>
#include <CGAL/IO/pixmaps/optimal_convex.xpm>
#include <CGAL/IO/pixmaps/points.xpm>

#include <qiconset.h>


Layers_toolbar::Layers_toolbar(CGAL::Qt_widget *w, QMainWindow *mw,
			       Cgal_Polygon *p) : QToolBar(mw, "LT"),
     nr_of_buttons(0)
  {
    showP = new Qt_layer_show_polygon<Cgal_Polygon>(*p);
    showGA = new Qt_layer_show_greene_approx<Cgal_Polygon>(*p);
    showYM = new Qt_layer_show_ymonotone<Cgal_Polygon>(*p);
    showOC = new Qt_layer_show_optimal_convex<Cgal_Polygon>(*p);
    showPP = new Qt_layer_show_polygon_points<Cgal_Polygon>(*p);

    //set the widget
    widget = w;
    window = mw;

    widget->attach(showP);
    widget->attach(showGA);
    widget->attach(showYM);
    widget->attach(showOC);
    widget->attach(showPP);

    showGA->deactivate();
    showYM->deactivate();
    showOC->deactivate();

    QIconSet set0(QPixmap( (const char**)show_polygon_small_xpm ),
                  QPixmap( (const char**)show_polygon_xpm ));
    QIconSet set1(QPixmap( (const char**)greene_approx_small_xpm ),
                  QPixmap( (const char**)greene_approx_xpm ));
    QIconSet set2(QPixmap( (const char**)ymonotone_small_xpm ),
                  QPixmap( (const char**)ymonotone_xpm ));
    QIconSet set3(QPixmap( (const char**)optimal_convex_small_xpm ),
                  QPixmap( (const char**)optimal_convex_xpm ));
    QIconSet set4(QPixmap( (const char**)points_small_xpm ),
                  QPixmap( (const char**)points_xpm ));

    but[0] = new QToolButton(this, "show_polygon");
    but[0]->setIconSet(set0);
    but[0]->setTextLabel("Show Simple Polygon");
    but[1] = new QToolButton(this, "greene_approx");
    but[1]->setIconSet(set1);
    but[1]->setTextLabel("Show Greene Approximation");
    but[2] = new QToolButton(this, "ymonotone");
    but[2]->setIconSet(set2);
    but[2]->setTextLabel("Show Y Monotone Partition");
    but[3] = new QToolButton(this, "optimal_convex");
    but[3]->setIconSet(set3);
    but[3]->setTextLabel("Show Optimal Convex Partition");
    but[4] = new QToolButton(this, "show_points");
    but[4]->setIconSet(set4);
    but[4]->setTextLabel("Show Polygon Vertices");

    nr_of_buttons = 5;
    button_group = new QButtonGroup(0, "nonexclusive");

    for(int i =0; i<nr_of_buttons; i++){
      but[i]->setToggleButton(true);
      but[i]->toggle();
      button_group->insert(but[i]);
    }
    but[1]->toggle();
    but[2]->toggle();
    but[3]->toggle();
    connect(but[0], SIGNAL(stateChanged(int)),
        showP, SLOT(stateChanged(int)));
    connect(but[1], SIGNAL(stateChanged(int)),
        showGA, SLOT(stateChanged(int)));
    connect(but[2], SIGNAL(stateChanged(int)),
        showYM, SLOT(stateChanged(int)));
    connect(but[3], SIGNAL(stateChanged(int)),
        showOC, SLOT(stateChanged(int)));
    connect(but[4], SIGNAL(stateChanged(int)),
        showPP, SLOT(stateChanged(int)));
    connect(button_group, SIGNAL(clicked(int)),
          widget, SLOT(redraw()));
  }

  Layers_toolbar::~Layers_toolbar()
  {
    delete showP;
    delete showGA;
    delete showYM;
    delete showOC;
    delete showPP;
    delete button_group;
  };

#include "partition_2_toolbar_layers.moc"

