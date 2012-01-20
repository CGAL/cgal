// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#include <CGAL/basic.h>


#include <CGAL/IO/Qt_widget.h>

#include "boolean_operations_2_toolbar.h"
#include "boolean_operations_2.h"

// icons
#include <CGAL/IO/pixmaps/arrow.xpm>
#include "icons/insert_circle.xpm"
#include "icons/insert_polygon.xpm"
#include "icons/locate.xpm"



#include <qiconset.h>


class MyWindow;
//mw used to be QMainWindow 
Tools_toolbar::Tools_toolbar(CGAL::Qt_widget *w,
                             MyWindow *mw) :
    QToolBar(mw, "NT")
  {


    w->attach(&getsimplebut);
    w->attach(&getcirclebut);
	 locatebut = new Qt_widget_locate_layer(mw);   
    //w->attach(&locatebut);
    w->attach(locatebut);
    
    /*w->attach(delete_red_but);
    w->attach(delete_blue_but);*/
    //w->attach(&delete_polygon);

    getsimplebut.deactivate();
    getcirclebut.deactivate();
    locatebut->deactivate();
    /*delete_red_but->deactivate();
    delete_blue_but->deactivate();*/
    //delete_polygon.deactivate();
    //set the widget
    widget = w;

    QIconSet set0(QPixmap( (const char**)arrow_small_xpm ),
                  QPixmap( (const char**)arrow_xpm ));
    QIconSet set1(QPixmap( (const char**)insert_polygon_xpm ),
                  QPixmap( (const char**)insert_polygon_xpm ));
    QIconSet set2(QPixmap( (const char**)insert_circle_xpm ),
                  QPixmap( (const char**)insert_circle_xpm ));
    QIconSet set3(QPixmap( (const char**)locate_xpm ),
                  QPixmap( (const char**)locate_xpm ));




  but[0] = new QToolButton(this, "deactivate layer");
  but[0]->setIconSet(set0);
  but[0]->setTextLabel("Deactivate Layer");
  but[1] = new QToolButton(this, "polygontool");
  but[1]->setIconSet(set1);
  but[1]->setTextLabel("Insert Circluar Polygon");
  but[2] = new QToolButton(this, "circletool");
  but[2]->setIconSet(set2);
  but[2]->setTextLabel("Insert Circle");
  but[3] = new QToolButton(this, "locatetool");
  but[3]->setIconSet(set3);
  but[3]->setTextLabel("Locate Polygon");


  int nr_of_buttons = 4;
  button_group = new QButtonGroup(0, "My_group");
  for(int i = 0; i < nr_of_buttons; i++)
  {
    button_group->insert(but[i]);
    but[i]->setToggleButton(true);
  }
  button_group->setExclusive(true);


  connect(but[1], SIGNAL(stateChanged(int)),
        &getsimplebut, SLOT(stateChanged(int)));

  connect(but[2], SIGNAL(stateChanged(int)),
        &getcirclebut, SLOT(stateChanged(int)));

  connect(but[3], SIGNAL(stateChanged(int)),
        locatebut, SLOT(stateChanged(int)));
}




#include "boolean_operations_2_toolbar.moc"

