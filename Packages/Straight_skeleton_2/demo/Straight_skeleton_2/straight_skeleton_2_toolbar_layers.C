// Copyright (c) 2002  Max Planck Institut fuer Informatik (Germany).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Radu Ursu

#ifdef CGAL_USE_QT

#include "straight_skeleton_2_toolbar_layers.h"
#include "straight_skeleton_2_layers.h"

// icons
#include <CGAL/IO/pixmaps/greene_approx.xpm>
#include <CGAL/IO/pixmaps/show_polygon.xpm>

#include <qiconset.h>

Layers_toolbar::Layers_toolbar(CGAL::Qt_widget *w
                              ,QMainWindow *mw
			      ,PolygonalRegion const& pr
			      ,Ssds const& ss
			      ) : QToolBar(mw, "LT"),
     nr_of_buttons(0)
  {
    showP  = new Qt_layer_show_polygon<PolygonalRegion>(pr);
    showSS = new Qt_layer_show_skeleton<Ssds>(ss);

    //set the widget
    widget = w;
    window = mw;

    widget->attach(showP);
    widget->attach(showSS);



    QIconSet set0(QPixmap( (const char**)show_polygon_small_xpm ),
                  QPixmap( (const char**)show_polygon_xpm ));
    QIconSet set1(QPixmap( (const char**)greene_approx_small_xpm ),
                  QPixmap( (const char**)greene_approx_xpm ));

    but[0] = new QToolButton(this, "show_polygon");
    but[0]->setIconSet(set0);
    but[0]->setTextLabel("Show Simple Polygon");
    but[1] = new QToolButton(this, "straight skeleton");
    but[1]->setIconSet(set1);
    but[1]->setTextLabel("Show Skeleton");
    
    nr_of_buttons = 2;
    button_group = new QButtonGroup(0, "nonexclusive");
    
    for(int i =0; i<nr_of_buttons; i++){
      but[i]->setToggleButton(true);
      but[i]->toggle();
      button_group->insert(but[i]);
    }
    //but[1]->toggle();
    connect(but[0], SIGNAL(stateChanged(int)),
        showP, SLOT(stateChanged(int)));
    connect(but[1], SIGNAL(stateChanged(int)),
        showSS, SLOT(stateChanged(int)));
  }	

  Layers_toolbar::~Layers_toolbar()
  {
    delete showP;
    delete showSS;
    delete button_group;
  };


#include "straight_skeleton_2_toolbar_layers.moc"

#endif
