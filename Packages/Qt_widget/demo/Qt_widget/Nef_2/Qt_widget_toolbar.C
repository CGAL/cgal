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
// file          : Qt_widget_toolbar.C
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau
//
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================

#if defined CGAL_USE_QT && defined CGAL_USE_GMP

#include <CGAL/IO/Qt_widget.h>
#include "Qt_widget_toolbar.h"

// icons
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/line.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>
#include <CGAL/IO/pixmaps/polygon.xpm>

#include <qiconset.h>

Tools_toolbar::Tools_toolbar(CGAL::Qt_widget *w, 
			     QMainWindow *mw) : QToolBar(mw, "NT")
  {
    w->attach(&input_point); 
    input_point.deactivate();
    w->attach(&input_line); 
    input_line.deactivate();
    w->attach(&input_polygon); 
    input_polygon.deactivate();
    
    //set the widget
    widget = w;

    QIconSet set0(QPixmap( (const char**)arrow_small_xpm ),
                  QPixmap( (const char**)arrow_xpm ));
    QIconSet set1(QPixmap( (const char**)point_small_xpm ),
                  QPixmap( (const char**)point_xpm ));
    QIconSet set2(QPixmap( (const char**)line_small_xpm ),
                  QPixmap( (const char**)line_xpm ));
    QIconSet set3(QPixmap( (const char**)polygon_small_xpm ),
                  QPixmap( (const char**)polygon_xpm ));
	
  but[0] = new QToolButton(this, "Deactivate Layer");
  but[0]->setIconSet(set0);
  but[0]->setTextLabel("Deactivate Layer");
  but[1] = new QToolButton(this, "pointinput layer");
  but[1]->setIconSet(set1);
  but[1]->setTextLabel("Input Point");
  but[2] = new QToolButton(this, "lineinput layer");
  but[2]->setIconSet(set2);
  but[2]->setTextLabel("Input Line");
  but[3] = new QToolButton(this, "polygoninput layer");
  but[3]->setIconSet(set3);
  but[3]->setTextLabel("Input Simple Polygon");
  
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

#include "Qt_widget_toolbar.moc"

#endif
