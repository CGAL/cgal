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


#include <CGAL/IO/Qt_widget.h>
#include "Qt_widget_toolbar.h"

// icons
#include <CGAL/IO/pixmaps/movepoint.xpm>
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>


Tools_toolbar::Tools_toolbar(CGAL::Qt_widget *w,
			     QMainWindow *mw, std::list<Point_2> *l1) :
  QToolBar(mw, "NT")
  {
    w->attach(&move_deletebut);
    w->attach(&pointbut);
    move_deletebut.deactivate();
    pointbut.deactivate();
    move_deletebut.pass_the_structure(l1);
    //set the widget
    widget = w;

    QIconSet set0(QPixmap( (const char**)arrow_small_xpm ),
                  QPixmap( (const char**)arrow_xpm ));
    QIconSet set1(QPixmap( (const char**)point_small_xpm ),
                  QPixmap( (const char**)point_xpm ));
    QIconSet set2(QPixmap( (const char**)movepoint_small_xpm ),
                  QPixmap( (const char**)movepoint_xpm ));

  but[0] = new QToolButton(this, "deactivate layer");
  but[0]->setIconSet(set0);
  but[0]->setTextLabel("Deactivate Layer");
  but[1] = new QToolButton(this, "pointtool");
  but[1]->setIconSet(set1);
  but[1]->setTextLabel("Input Point");
  but[2] = new QToolButton(this, "move/delete tool");
  but[2]->setIconSet(set2);
  but[2]->setTextLabel("Move/Delete Point");

  nr_of_buttons = 3;
  button_group = new QButtonGroup(0, "My_group");
  for(int i = 0; i<nr_of_buttons; i++) {
    button_group->insert(but[i]);
    but[i]->setToggleButton(true);
  }
  button_group->setExclusive(true);

  connect(but[1], SIGNAL(stateChanged(int)),
        &pointbut, SLOT(stateChanged(int)));
  connect(but[2], SIGNAL(stateChanged(int)),
        &move_deletebut, SLOT(stateChanged(int)));
}

#include "Qt_widget_toolbar.moc"

