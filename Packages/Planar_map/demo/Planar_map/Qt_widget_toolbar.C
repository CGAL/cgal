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
// file          : demo/Qt_widget/Planar_map/Qt_widget_toolbar.C
// package       : Planar_map
// author(s)     : Efi Fogel
// release       : 
// release_date  : 
//
// coordinator   : Efi Fogel
//
// ============================================================================

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/pixmaps/movepoint.xpm>
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/line.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>

#include <qiconset.h>

#include "Qt_widget_toolbar.h"

Tools_toolbar::Tools_toolbar(CGAL::Qt_widget * w, 
			     QMainWindow *mw, std::list<Curve> * l1, Planar_map * pm) :
  QToolBar(mw, "NT")
{
  snapping_layer.pass_the_structure(l1, pm);
  segment_layer.pass_the_structure(l1, pm);
  w->attach(&point_layer);
  w->attach(&segment_layer);
  w->attach(&snapping_layer);

  //    move_deletebut.deactivate();
  point_layer.deactivate();
  segment_layer.deactivate();
  snapping_layer.deactivate();

  //    move_deletebut.pass_the_structure(l1);
  //set the widget
  widget = w;

  QIconSet set0(QPixmap( (const char**)arrow_small_xpm ),
                QPixmap( (const char**)arrow_xpm ));
  QIconSet set1(QPixmap( (const char**)point_small_xpm ),
                QPixmap( (const char**)point_xpm ));
  QIconSet set2(QPixmap( (const char**)line_small_xpm ),
                QPixmap( (const char**)line_xpm ));
  QIconSet set3(QPixmap( (const char**)movepoint_small_xpm ),
                QPixmap( (const char**)movepoint_xpm ));
    

  but[0] = new QToolButton(this, "deactivate layer");
  but[0]->setIconSet(set0);
  but[0]->setTextLabel("Deactivate Layer");
  but[1] = new QToolButton(this, "pointinput layer");
  but[1]->setIconSet(set1);
  but[1]->setTextLabel("Input point");
  but[2] = new QToolButton(this, "segment layer");
  but[2]->setIconSet(set2);
  but[2]->setTextLabel("Segment layer");
  but[3] = new QToolButton(this, "snapping layer");
  but[3]->setIconSet(set3);
  but[3]->setTextLabel("Snapping layer");

  nr_of_buttons = 4;
  button_group = new QButtonGroup(0, "My_group");
  for (int i = 0; i<nr_of_buttons; i++) {
    button_group->insert(but[i]);
    but[i]->setToggleButton(true);
  }
  button_group->setExclusive(true);
  
  connect(but[1], SIGNAL(stateChanged(int)),
          &point_layer, SLOT(stateChanged(int)));
  connect(but[2], SIGNAL(stateChanged(int)),
          &segment_layer, SLOT(stateChanged(int)));
  connect(but[3], SIGNAL(stateChanged(int)),
          &snapping_layer, SLOT(stateChanged(int)));

};

#include "Qt_widget_toolbar.moc"

#endif
