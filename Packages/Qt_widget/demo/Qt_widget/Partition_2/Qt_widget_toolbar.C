// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium

// This software and related documentation are part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation are provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind. 
//
// Every use of CGAL requires a license. 
//
// Academic research and teaching license
// - For academic research and teaching purposes, permission to use and copy
//   the software and its documentation is hereby granted free of charge,
//   provided that it is not a component of a commercial product, and this
//   notice appears in all copies of the software and related documentation. 
//
// Commercial licenses
// - Please check the CGAL web site http://www.cgal.org/index2.html for 
//   availability.
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// file          : src/Qt_widget_toolbar.C
// package       : Qt_widget
// author(s)     : Ursu Radu
// release       : CGAL-2.4
// release_date  : 2002, May 16
//
// coordinator   : Laurent Rineau
//
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================

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
