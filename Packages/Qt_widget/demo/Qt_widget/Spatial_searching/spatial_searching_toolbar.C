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
// file          : spatial_searching_toolbar.C
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
#include "spatial_searching_toolbar.h"

// icons
#include <CGAL/IO/pixmaps/movepoint.xpm>
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>
#include <CGAL/IO/pixmaps/circle.xpm>
#include <CGAL/IO/pixmaps/iso_rectangle.xpm>


#include <qiconset.h>

Tools_toolbar::Tools_toolbar(CGAL::Qt_widget *w, 
			     QMainWindow *mw, std::vector<Point> *l1) :
  QToolBar(mw, "NT")
  {
    w->attach(&edit_layer);
    w->attach(&point_layer);
    w->attach(&iso_r_layer);
    w->attach(&circle_layer);
    edit_layer.deactivate();
    point_layer.deactivate();
    iso_r_layer.deactivate();
    circle_layer.deactivate();
    edit_layer.pass_the_structure(l1);

    //set the widget
    widget = w;

    QIconSet set0(QPixmap( (const char**)arrow_small_xpm ),
                  QPixmap( (const char**)arrow_xpm ));
    QIconSet set1(QPixmap( (const char**)point_small_xpm ),
                  QPixmap( (const char**)point_xpm ));
    QIconSet set2(QPixmap( (const char**)iso_rectangle_small_xpm ),
                  QPixmap( (const char**)iso_rectangle_xpm ));
    QIconSet set3(QPixmap( (const char**)circle_small_xpm ),
                  QPixmap( (const char**)circle_xpm ));
    QIconSet set4(QPixmap( (const char**)movepoint_small_xpm ),
                  QPixmap( (const char**)movepoint_xpm ));


  but[0] = new QToolButton(this, "deactivate layer");
  but[0]->setIconSet(set0);
  but[0]->setTextLabel("Deactivate Layer");
  but[1] = new QToolButton(this, "pointtool");
  but[1]->setIconSet(set1);
  but[1]->setTextLabel("Input Point");
  but[2] = new QToolButton(this, "iso_rectangle");
  but[2]->setIconSet(set2);
  but[2]->setTextLabel("Input Iso_rectangle");
  but[3] = new QToolButton(this, "iso_rectangle");
  but[3]->setIconSet(set3);
  but[3]->setTextLabel("Input Circle");
  but[4] = new QToolButton(this, "move/delete tool");
  but[4]->setIconSet(set4);
  but[4]->setTextLabel("Move/Delete Point");
  
  nr_of_buttons = 5;
  button_group = new QButtonGroup(0, "My_group");
  for(int i = 0; i<nr_of_buttons; i++) {
    button_group->insert(but[i]);
    but[i]->setToggleButton(true);
  }
  button_group->setExclusive(true);
  
  connect(but[1], SIGNAL(stateChanged(int)),
        &point_layer, SLOT(stateChanged(int)));
  connect(but[2], SIGNAL(stateChanged(int)),
        &iso_r_layer, SLOT(stateChanged(int)));
  connect(but[3], SIGNAL(stateChanged(int)),
        &circle_layer, SLOT(stateChanged(int)));
  connect(but[4], SIGNAL(stateChanged(int)),
        &edit_layer, SLOT(stateChanged(int)));
};

#include "spatial_searching_toolbar.moc"

#endif
