// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : demo/Qt_widget/Max_k-gon/Qt_widget_toolbar.C
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_widget.h>
#include "Qt_widget_toolbar.h"

// icons
#include <CGAL/IO/pixmaps/movepoint.xpm>
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>


namespace CGAL {
  Tools_toolbar::Tools_toolbar(Qt_widget *w, 
			  QMainWindow *mw, std::list<Point> *l1)
  {
    w->attach(&move_deletebut);
    w->attach(&pointbut);
    move_deletebut.deactivate();
    pointbut.deactivate();
    move_deletebut.pass_the_structure(l1);
    //set the widget
    widget = w;

#if QT_VERSION < 300
  // for Qt 2.3 and before
  maintoolbar = new QToolBar("tools", mw, 
      QMainWindow::Top, TRUE, "Tools");
#else
  // from Qt 3.0
  maintoolbar = new QToolBar(mw, "Tools");
  mw->addDockWindow (maintoolbar, "tools", DockTop, TRUE );
#endif
		
  but[0] = new QToolButton(maintoolbar, "notool");
  but[0]->setPixmap(QPixmap( (const char**)arrow_xpm ));
  but[1] = new QToolButton(maintoolbar, "pointtool");
  but[1]->setPixmap(QPixmap( (const char**)point_xpm ));
  but[2] = new QToolButton(maintoolbar, "move/delete tool");
  but[2]->setPixmap(QPixmap( (const char**)movepoint_xpm ));
  
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
};
        
}//end namespace CGAL

#include "Qt_widget_toolbar.moc"

#endif
