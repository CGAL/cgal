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
// file          : demo/Qt_widget/Nef_2/Qt_widget_toolbar.C
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : CGAL-2.4
// release_date  : 2002, July 8
//
// coordinator   : Laurent Rineau
//
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_widget.h>
#include "Qt_widget_toolbar.h"

// icons
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/line.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>
#include <CGAL/IO/pixmaps/polygon.xpm>


namespace CGAL {
  Tools_toolbar::Tools_toolbar(Qt_widget *w, 
			  QMainWindow *mw)
  {
    w->attach(&input_point); 
    input_point.deactivate();
    w->attach(&input_line); 
    input_line.deactivate();
    w->attach(&input_polygon); 
    input_polygon.deactivate();
    
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
  but[2] = new QToolButton(maintoolbar, "linetool");
  but[2]->setPixmap(QPixmap( (const char**)line_xpm ));
  but[3] = new QToolButton(maintoolbar, "polygontool");
  but[3]->setPixmap(QPixmap( (const char**)polygon_xpm ));
  
  nr_of_buttons = 4;
  button_group = new QButtonGroup(0, "My_group");
  for(int i = 0; i<nr_of_buttons; i++) {
    button_group->insert(but[i]);
    but[i]->setToggleButton(true);
  }
  button_group->setExclusive(true);
  
  connect(but[1], SIGNAL(stateChanged(int)),
        &input_point, SLOT(stateChanged(int)));
  connect(but[2], SIGNAL(stateChanged(int)),
        &input_line, SLOT(stateChanged(int)));
  connect(but[3], SIGNAL(stateChanged(int)),
        &input_polygon, SLOT(stateChanged(int)));
};
        
}//end namespace CGAL

#include "Qt_widget_toolbar.moc"

#endif
