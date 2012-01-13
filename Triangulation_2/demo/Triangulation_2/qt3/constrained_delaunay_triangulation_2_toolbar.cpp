// Copyright (c) 1997-2002  INRIA Sophia-Antipolis (France).
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


#include <CGAL/IO/Qt_widget.h>
#include "constrained_delaunay_triangulation_2_toolbar.h"

// icons
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>
#include <CGAL/IO/pixmaps/line.xpm>
#include <CGAL/IO/pixmaps/polygon.xpm>

#include <qiconset.h>

  Tools_toolbar::Tools_toolbar(CGAL::Qt_widget *w, QMainWindow *mw,
			       CDT* ) : QToolBar(mw, "NT")
  {
    //when it is created, the toolbar has 0 buttons
    nr_of_buttons = 0;
    //set the widget
    widget = w;
    widget->attach(&segmentbut);
    widget->attach(&pointbut);
    widget->attach(&polygonbut);
    pointbut.deactivate();
    segmentbut.deactivate();
    polygonbut.deactivate();

    QIconSet set0(QPixmap( (const char**)arrow_small_xpm ),
                  QPixmap( (const char**)arrow_xpm ));
    QIconSet set1(QPixmap( (const char**)point_small_xpm ),
                  QPixmap( (const char**)point_xpm ));
    QIconSet set2(QPixmap( (const char**)line_small_xpm ),
                  QPixmap( (const char**)line_xpm ));
    QIconSet set3(QPixmap( (const char**)polygon_small_xpm ),
                  QPixmap( (const char**)polygon_xpm ));
    but[0] = new QToolButton(this, "deactivate layer");
    but[0]->setIconSet(set0);
    but[0]->setTextLabel("Deactivate Layer");
    but[1] = new QToolButton(this, "pointinput layer");
    but[1]->setIconSet(set1);
    but[1]->setTextLabel("Input Point");
    but[2] = new QToolButton(this, "lineinput layer");
    but[2]->setIconSet(set2);
    but[2]->setTextLabel("Input Segment");
    but[3] = new QToolButton(this, "polygoninput layer");
    but[3]->setIconSet(set3);
    but[3]->setTextLabel("Input Polygon");

  nr_of_buttons = 4;

  button_group = new QButtonGroup(0, "exclusive");
  for (int i=0; i<nr_of_buttons; i++) {
    button_group->insert(but[i]);
    but[i]->setToggleButton(true);
  }
  button_group->setExclusive(true);
  connect(but[1], SIGNAL(stateChanged(int)),
        &pointbut, SLOT(stateChanged(int)));
  connect(but[2], SIGNAL(stateChanged(int)),
        &segmentbut, SLOT(stateChanged(int)));
  connect(but[3], SIGNAL(stateChanged(int)),
        &polygonbut, SLOT(stateChanged(int)));
  but[2]->toggle();

};

#include "constrained_delaunay_triangulation_2_toolbar.moc"

