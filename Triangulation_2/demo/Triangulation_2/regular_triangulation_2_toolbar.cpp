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


#include <CGAL/IO/Qt_widget.h>
#include "regular_triangulation_2_toolbar.h"

// icons
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/movepoint.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>
#include <CGAL/IO/pixmaps/circle.xpm>

#include <qiconset.h>

  Tools_toolbar::Tools_toolbar(CGAL::Qt_widget *w, QMainWindow *mw,
			       Regular_triangulation *t) :
    QToolBar(mw, "NT")
  {
    //when it is created, the toolbar has 0 buttons
    nr_of_buttons = 0;
    //set the widget
    widget = w;
    widget->attach(&input_circle_layer);
    widget->attach(&input_point_layer);
    widget->attach(&edit_vertex_layer);
    input_point_layer.deactivate();
    input_circle_layer.deactivate();
    edit_vertex_layer.set_triangulation(t);
    edit_vertex_layer.deactivate();

    QIconSet set0(QPixmap( (const char**)arrow_small_xpm ),
                  QPixmap( (const char**)arrow_xpm ));
    QIconSet set1(QPixmap( (const char**)point_small_xpm ),
                  QPixmap( (const char**)point_xpm ));
    QIconSet set2(QPixmap( (const char**)circle_small_xpm ),
                  QPixmap( (const char**)circle_xpm ));
    QIconSet set3(QPixmap( (const char**)movepoint_small_xpm ),
                  QPixmap( (const char**)movepoint_xpm ));

    but[0] = new QToolButton(this, "deactivate layer");
    but[0]->setIconSet(set0);
    but[1] = new QToolButton(this, "pointinput layer");
    but[1]->setIconSet(set1);
    but[1]->setTextLabel("Input Point");
    but[2] = new QToolButton(this, "circleinput layer");
    but[2]->setIconSet(set2);
    but[2]->setTextLabel("Input WeightedPoint");
    but[3] = new QToolButton(this, "movedelete layer");
    but[3]->setIconSet(set3);
    but[3]->setTextLabel("Edit Vertex");

  nr_of_buttons = 4;

  button_group = new QButtonGroup(0, "exclusive");
  for (int i=0; i<nr_of_buttons; i++) {
    button_group->insert(but[i]);
    but[i]->setToggleButton(true);
  }
  button_group->setExclusive(true);
  connect(but[1], SIGNAL(stateChanged(int)),
        &input_point_layer, SLOT(stateChanged(int)));
  connect(but[2], SIGNAL(stateChanged(int)),
        &input_circle_layer, SLOT(stateChanged(int)));
  connect(but[3], SIGNAL(stateChanged(int)),
        &edit_vertex_layer, SLOT(stateChanged(int)));
  but[2]->toggle();

};



#include "regular_triangulation_2_toolbar.moc"

