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
// file          : src/Qt_widget_toolbar.C
// package       : Qt_widget
// author(s)     : Ursu Radu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifdef CGAL_USE_QT

#include "Qt_widget_toolbar.h"
//#include <CGAL/IO/Qt_widget.h>

// icons

#include <CGAL/IO/pixmaps/arrow.xpm>
#include <CGAL/IO/pixmaps/polygon.xpm>


namespace CGAL {
  Tools_toolbar::Tools_toolbar(Qt_widget *w, QMainWindow *mw)
  {
    //when it is created, the toolbar has 0 buttons
    nr_of_buttons = 0;
    //set the widget
    widget = w;
    widget->attach(&getsimplebut);

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
  
  but[1] = new QToolButton(maintoolbar, "spolygon");
  but[1]->setPixmap(QPixmap( (const char**)polygon_xpm ));
  but[1]->setTextLabel("Input Simple Polygon");
  
  button_group = new QButtonGroup(0, "exclusive_group");
  button_group->insert(but[0]);
  button_group->insert(but[1]);
  button_group->setExclusive(true);
  
  but[1]->setToggleButton(TRUE);
  
  connect(but[1], SIGNAL(stateChanged(int)),
        &getsimplebut, SLOT(stateChanged(int)));
  nr_of_buttons = 2;
  };

	

}//end namespace

#include "Qt_widget_toolbar.moc"

#endif
