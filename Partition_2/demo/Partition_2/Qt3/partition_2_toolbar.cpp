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


#include "partition_2_toolbar.h"

// icons

#include <CGAL/IO/pixmaps/arrow.xpm>
#include <CGAL/IO/pixmaps/polygon.xpm>

#include <qiconset.h>

Tools_toolbar::Tools_toolbar(CGAL::Qt_widget *w, QMainWindow *mw) :
  QToolBar(mw, "NT")
{
    //when it is created, the toolbar has 0 buttons
    nr_of_buttons = 0;
    //set the widget
    widget = w;
    widget->attach(&getsimplebut);
    getsimplebut.deactivate();

    QIconSet set0(QPixmap( (const char**)arrow_small_xpm ),
                  QPixmap( (const char**)arrow_xpm ));
    QIconSet set1(QPixmap( (const char**)polygon_small_xpm ),
                  QPixmap( (const char**)polygon_xpm ));

  but[0] = new QToolButton(this, "deactivate layer");
  but[0]->setIconSet(set0);
  but[0]->setTextLabel("Deactivate Layer");

  but[1] = new QToolButton(this, "spolygon");
  but[1]->setIconSet(set1);
  but[1]->setTextLabel("Input Simple Polygon");

  button_group = new QButtonGroup(0, "exclusive_group");
  button_group->insert(but[0]);
  button_group->insert(but[1]);
  button_group->setExclusive(true);

  but[0]->setToggleButton(true);
  but[1]->setToggleButton(true);

  connect(but[1], SIGNAL(stateChanged(int)),
        &getsimplebut, SLOT(stateChanged(int)));
  nr_of_buttons = 2;
  };

#include "partition_2_toolbar.moc"

