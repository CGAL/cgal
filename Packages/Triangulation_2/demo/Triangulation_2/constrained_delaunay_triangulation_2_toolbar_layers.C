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
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// file          : constrained_delaunay_triangulation_2_toolbar_layers.C
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


#ifdef CGAL_USE_QT

#include "constrained_delaunay_triangulation_2_toolbar_layers.h"

// icons
#include <CGAL/IO/pixmaps/points.xpm>
#include <CGAL/IO/pixmaps/constrained.xpm>
#include <CGAL/IO/pixmaps/triangulation.xpm>

#include <qiconset.h>

  Layers_toolbar::Layers_toolbar(CGAL::Qt_widget *w, QMainWindow *mw, 
                                 CDT *t) : QToolBar(mw, "LT"), nr_of_buttons(0)
  {
    showT   = new Qt_layer_show_triangulation< CDT >(*t);
    showP   = new Qt_layer_show_points< CDT >(*t);
    showC   = new Qt_layer_show_constraints<CDT>(*t);

    //set the widget
    widget = w;
    window = mw;

    widget->attach(showT);    
    widget->attach(showC);
    widget->attach(showP);
    
    QIconSet set0(QPixmap( (const char**)triangulation_small_xpm ),
                  QPixmap( (const char**)triangulation_xpm ));
    QIconSet set1(QPixmap( (const char**)constrained_small_xpm ),
                  QPixmap( (const char**)constrained_xpm ));
    QIconSet set2(QPixmap( (const char**)points_small_xpm ),
                  QPixmap( (const char**)points_xpm ));
		
    but[0] = new QToolButton(this, "triangulation");
    but[0]->setIconSet(set0);
    but[0]->setTextLabel("Show Triangulation");
    but[1] = new QToolButton(this, "constraints");
    but[1]->setIconSet(set1);
    but[1]->setTextLabel("Show Constraints");
    but[2] = new QToolButton(this, "vertices");
    but[2]->setIconSet(set2);
    but[2]->setTextLabel("Show Vertices");
		

    nr_of_buttons = 3;
	  button_group = new QButtonGroup(0, "nonexclusive");
    for(int i =0; i<nr_of_buttons; i++)
    {
      but[i]->setToggleButton(TRUE);
      button_group->insert(but[i]);
      but[i]->toggle();
    }    
    connect(button_group, SIGNAL(clicked(int)),
          widget, SLOT(redraw()));
    
    connect(but[0], SIGNAL(stateChanged(int)),
        showT, SLOT(stateChanged(int)));
    connect(but[1], SIGNAL(stateChanged(int)),
        showC, SLOT(stateChanged(int)));
    connect(but[2], SIGNAL(stateChanged(int)),
        showP, SLOT(stateChanged(int)));
  }
  

#include "constrained_delaunay_triangulation_2_toolbar_layers.moc"

#endif
