// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
#include <CGAL/IO/Qt_widget_standard_toolbar.h>

// icons
#include <CGAL/IO/pixmaps/zoom_in_rect.xpm>
#include <CGAL/IO/pixmaps/zoom_out.xpm>
#include <CGAL/IO/pixmaps/zoom_in.xpm>
#include <CGAL/IO/pixmaps/focus.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>
#include <CGAL/IO/pixmaps/back.xpm>
#include <CGAL/IO/pixmaps/forward.xpm>
#include <CGAL/IO/pixmaps/mouse_coord.xpm>

#include <qbuttongroup.h>
#include <qiconset.h>
#include <qmainwindow.h>
#include <qtoolbutton.h>
#include <qstring.h>
#include <qvariant.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_history.h>
#include <CGAL/IO/Qt_widget_focus.h>
#include <CGAL/IO/Qt_widget_zoomrect.h>
#include <CGAL/IO/Qt_widget_handtool.h>
#include <CGAL/IO/Qt_widget_show_mouse_coordinates.h>

namespace CGAL {
  Qt_widget_standard_toolbar::
  Qt_widget_standard_toolbar(Qt_widget *w, QMainWindow *parent,
			     const char* name) :
    QToolBar(parent, name),
    widget(w)
  {
    setLabel("Qt_widget standard toolbar");
    fill_toolbar(parent);
  }
  
  Qt_widget_standard_toolbar::
  Qt_widget_standard_toolbar(Qt_widget *w, QMainWindow *mw,
			     QWidget* parent,
			     bool newLine,
			     const char* name) :
    QToolBar("Qt_widget standard toolbar", mw, parent, newLine, name),
    widget(w)
  {
    fill_toolbar(mw);
  }
  
  void Qt_widget_standard_toolbar::fill_toolbar(QMainWindow *mw)
  {
    Qt_widget_focus* focuslayer = 
      new Qt_widget_focus(this, "focuslayer");
    Qt_widget_zoomrect* zoomrectlayer = 
      new Qt_widget_zoomrect(this, "zoomrectlayer");
    Qt_widget_handtool* handtoollayer = 
      new Qt_widget_handtool(this, "handtoollayer");
    Qt_widget_show_mouse_coordinates* showcoordlayer = 0; // created below

    widget->attach_standard(focuslayer);
    widget->attach_standard(zoomrectlayer);
    widget->attach_standard(handtoollayer);
    focuslayer->deactivate();
    zoomrectlayer->deactivate();
    handtoollayer->deactivate();

    if (mw)
      {
	mw->statusBar();
	showcoordlayer = 
	  new Qt_widget_show_mouse_coordinates(*mw, this,
					       "showcoordlayer");
	widget->attach_standard(showcoordlayer);
	showcoordlayer->does_eat_events = false;
	widget->setMouseTracking(true);
      }

    const QPixmap p1(arrow_small_xpm), p2(arrow_xpm);
    const QIconSet arrow_pixmap(p1, p2);

    const QPixmap p3(back_small_xpm), p4(back_xpm);
    const QIconSet back_pixmap(p3, p4);

    const QPixmap p5(forward_small_xpm), p6(forward_xpm);
    const QIconSet forward_pixmap(p5, p6);

    QToolButton* backBt = new QToolButton(this, "History Back");
    backBt->setIconSet(back_pixmap);
    backBt->setTextLabel("History Back");

    QToolButton* forwardBt = new QToolButton(this, "History Forward");
    forwardBt->setIconSet(forward_pixmap);
    forwardBt->setTextLabel("History Forward");
  
    addSeparator();

    QToolButton* zoominBt =
      new QToolButton(QPixmap(zoomin_xpm ),
			       "Zoom in", 
			       0, 
			       this, 
			       SLOT(zoomin()), 
			       this, 
			       "Zoom in");
    zoominBt->setTextLabel("Scaling factor X2");

    QToolButton* zoomoutBt = 
      new QToolButton(QPixmap(zoomout_xpm ),
				"Zoom out", 
				0, 
				this, 
				SLOT(zoomout()), 
				this, 
				"Zoom out");
    zoomoutBt->setTextLabel("Scaling factor 1/2");

    addSeparator();

    nolayerBt = new QToolButton(this, "nolayer");
    nolayerBt->setIconSet(arrow_pixmap);
    nolayerBt->setTextLabel("Deactivate Standard Layer");

    QToolButton* zoomrectBt = new QToolButton(this, "focus on region");
    zoomrectBt->setPixmap(QPixmap(zoomin_rect_xpm ));
    zoomrectBt->setTextLabel("Focus on region");

    QToolButton* focusBt = new QToolButton(this, "focus");
    focusBt->setPixmap(QPixmap(focus_xpm ));
    focusBt->setTextLabel("Focus on point");

    QToolButton* handtoolBt = new QToolButton(this, "handtool");
    handtoolBt->setPixmap(QPixmap(hand_xpm ));
    handtoolBt->setTextLabel("Pan tool");

    addSeparator();

    QToolButton* showcoordBt = new QToolButton(this, "mouse");
    showcoordBt->setPixmap(QPixmap(mouse_coord_xpm) );
    showcoordBt->setTextLabel("Mouse Coordinates");

    button_group = new QButtonGroup(0, "My_group");
    // this button has no parent and is destroyed manually in the
    // destructor

    // below is the list of buttons in the group
    QToolButton* const button_group_list[] = { nolayerBt,
					       zoomrectBt,
					       focusBt,
					       handtoolBt };
    for(int i=0; i<4; ++i)
      {
	button_group_list[i]->setToggleButton(true);
	button_group->insert(button_group_list[i]);
      }
    button_group->setExclusive(true);
    connect(button_group, SIGNAL(clicked(int)),
	    this, SLOT(group_clicked(int)));

    nolayerBt->setOn(true);

    showcoordBt->setToggleButton(true);
    showcoordBt->toggle();

    connect(zoomrectBt, SIGNAL(stateChanged(int)),
	    zoomrectlayer, SLOT(stateChanged(int)));
    connect(focusBt, SIGNAL(stateChanged(int)),
	    focuslayer, SLOT(stateChanged(int)));
    connect(handtoolBt, SIGNAL(stateChanged(int)),
	    handtoollayer, SLOT(stateChanged(int)));
    connect(showcoordBt, SIGNAL(stateChanged(int)),
	    showcoordlayer, SLOT(stateChanged(int)));


    // history setting
    history = new Qt_widget_history(widget, "standard history");
    connect(backBt, SIGNAL(clicked()),
	    history, SLOT(backward()));
    connect(forwardBt, SIGNAL(clicked()),
	    history, SLOT(forward()));

    connect(history, SIGNAL(backwardAvaillable(bool)),
            backBt, SLOT(setEnabled(bool)));
    connect(history, SIGNAL(forwardAvaillable(bool)),
            forwardBt, SLOT(setEnabled(bool)));
    history->clear();
  }

  void Qt_widget_standard_toolbar::group_clicked(int i)
  {
    static int id = 0; 
    // This id is here to keep track of the button from the group that
    // is on (if all toolbuttons are down, nolayer is on. At the
    // beginning, it is set to 0, because
    // button_group.id(nolayerBt)==0.

    if( i == id )
      {
	if( i == 0) return;
	// nolayer is on and cannot be set off like that.

	QToolButton* tBt = 
	  static_cast<QToolButton*>(button_group->find(i));
	if( tBt != 0)
	  tBt->setOn(false);

	nolayerBt->setOn(true);
	id = 0;
      }
    else
      id = i;
  }

  void Qt_widget_standard_toolbar::zoomin()
  {
    widget->zoom(2); //
  }
  void Qt_widget_standard_toolbar::zoomout()
  {
    widget->zoom(0.5);
  }
  void Qt_widget_standard_toolbar::back()
  {
    history->backward();
  }
  void Qt_widget_standard_toolbar::forward()
  {
    history->forward();
  }
  void Qt_widget_standard_toolbar::clear_history()
  {
    history->clear();
  }

}//end namespace
#include "Qt_widget_standard_toolbar.moc"
