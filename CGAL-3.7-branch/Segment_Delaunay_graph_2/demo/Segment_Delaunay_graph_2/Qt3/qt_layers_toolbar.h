// Copyright (c) 2003,2005  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>

#ifndef QT_LAYERS_TOOLBAR_H
#define QT_LAYERS_TOOLBAR_H

#include <CGAL/IO/Qt_widget.h>

#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qstatusbar.h>
#include <qstring.h>
#include <qwhatsthis.h>

#include "qt_layers.h"

// icons
#include <CGAL/IO/pixmaps/points.xpm>
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/line.xpm>
#include <CGAL/IO/pixmaps/polygon.xpm>
#include <CGAL/IO/pixmaps/voronoi.xpm>
#include <CGAL/IO/pixmaps/triangulation.xpm>
#include <CGAL/IO/pixmaps/mouse_coord.xpm>
#include <CGAL/IO/pixmaps/notool.xpm>
#include <CGAL/IO/pixmaps/holddown.xpm>
//#include "remove.xpm"

typedef enum { SDG_POINT, SDG_SEGMENT, SDG_POLYGON } Input_mode;

class Layers_toolbar : public QToolBar
{
  Q_OBJECT
public:
  Layers_toolbar(CGAL::Qt_widget *widget, SDG_2& sdg,
		 const QString& label, QMainWindow* mainWindow,
		 QWidget* parent, bool newLine = FALSE,
		 const char* name = 0, WFlags f = 0, bool is_pvd = false)
    : QToolBar(label, mainWindow, parent, newLine, name, f),
      nr_of_buttons(0), input_mode(SDG_SEGMENT)
  {
    showVD = new Voronoi_diagram_layer<SDG_2>(sdg);
    showSI = new Sites_layer<SDG_2>(sdg);
    showSK = new Skeleton_layer<SDG_2>(sdg);

    // set the widget
    this->widget = widget;
    window = mainWindow;
    window->statusBar();

    widget->attach(showVD);
    widget->attach(showSK);
    widget->attach(showSI);


    but[0] = new QToolButton(QPixmap( (const char**)points_xpm ),
			     "Show sites",
			     0,
			     this,
			     SLOT(show_sites()),
			     this,
			     "Show sites");

    but[1] = new QToolButton(QPixmap( (const char**)voronoi_xpm ),
			     "Show Voronoi diagram",
			     0,
			     this,
			     SLOT(show_voronoi()),
			     this,
			     "Show Voronoi diagram");

    but[2] = new QToolButton(QPixmap( (const char**)triangulation_xpm ),
			     "Show skeleton",
			     0,
			     this,
			     SLOT(show_skeleton()),
			     this,
			     "Show skeleton");

    but[3] = new QToolButton(QPixmap( (const char**)point_xpm ),
			     "Input points",
			     0,
			     this,
			     SLOT(input_points_mode()),
			     this,
			     "Input points");


    but[4] = new QToolButton(QPixmap( (const char**)line_xpm ),
			     "Input segments",
			     0,
			     this,
			     SLOT(input_segments_mode()),
			     this,
			     "Input segments");
    but[5] = new QToolButton(QPixmap( (const char**)polygon_xpm ),
			     "Input polygons",
			     0,
			     this,
			     SLOT(input_polygons_mode()),
			     this,
			     "Input polygons");
#if 1
    but[6] = new QToolButton(QPixmap( (const char**)notool_xpm ),
			     "Remove site",
			     0,
			     this,
			     SLOT(remove_mode()),
			     this,
			     "Remove site");
#else
    but[6] = new QToolButton(QPixmap( (const char**)remove_xpm ),
			     "Remove site",
			     0,
			     this,
			     SLOT(remove_mode()),
			     this,
			     "Remove site");
#endif
    but[7] = new QToolButton(QPixmap( (const char**)holddown_xpm ),
			     "Snap mode",
			     0,
			     this,
			     SLOT(snap_mode()),
			     this,
			     "Snap mode");

    showSK->deactivate();

    nr_of_buttons = 8;
    for(int i = 0; i < nr_of_buttons; i++){
      but[i]->setToggleButton(TRUE);
    }

    but[0]->toggle();
    but[1]->toggle();
    //    but[2]->toggle();
    //    but[3]->toggle();
    if ( is_pvd ) {
      but[5]->toggle();
      input_mode = SDG_POLYGON;
    } else {
      but[4]->toggle();
      input_mode = SDG_SEGMENT;
    }
    //    but[6]->toggle();
    //    but[7]->toggle();
  }

  ~Layers_toolbar() {
    delete showVD;
    delete showSI;
    delete showSK;
  }

  inline QToolBar* toolbar() { return this; };

signals:
  void new_object(CGAL::Object);
  void inputModeChanged(Input_mode);
  void insertModeChanged(bool);
  void snapModeChanged(bool);

private slots:
  void show_sites() {
    if ( but[0]->isOn() ) {
      showSI->activate();
    } else {
      showSI->deactivate();
    }
    widget->redraw();
  }

  void show_voronoi() {
    if ( !but[1]->isOn() ) {
      showVD->deactivate();
    } else {
      showVD->activate();
      if ( but[2]->isOn() ) {
	but[2]->toggle();
	showSK->deactivate();
      }
    }

    widget->redraw();
  }

  void show_skeleton() {
    if ( !but[2]->isOn() ) {
      showSK->deactivate();
    } else {
      showSK->activate();
      if ( but[1]->isOn() ) {
	but[1]->toggle();
	showVD->deactivate();
      }
    }

    widget->redraw();
  }

  void input_points_mode() {
    if ( !but[3]->isOn() ) {
      but[3]->toggle();
      return;
    }
    if ( input_mode == SDG_SEGMENT ) { but[4]->toggle(); }
    if ( input_mode == SDG_POLYGON ) { but[5]->toggle(); }

    input_mode = SDG_POINT;
    emit inputModeChanged( SDG_POINT );
  }

  void input_segments_mode() {
    if ( !but[4]->isOn() ) {
      but[4]->toggle();
      return;
    }
    if ( input_mode == SDG_POINT ) { but[3]->toggle(); }
    if ( input_mode == SDG_POLYGON ) { but[5]->toggle(); }

    input_mode = SDG_SEGMENT;
    emit inputModeChanged( SDG_SEGMENT );
  }

  void input_polygons_mode() {
    if ( !but[5]->isOn() ) {
      but[5]->toggle();
      return;
    }
    if ( input_mode == SDG_POINT ) { but[3]->toggle(); }
    if ( input_mode == SDG_SEGMENT ) { but[4]->toggle(); }

    input_mode = SDG_POLYGON;
    emit inputModeChanged( SDG_POLYGON );
  }

  void remove_mode() {
    emit insertModeChanged( but[6]->isOn() );
  }

  void snap_mode() {
    emit snapModeChanged( but[7]->isOn() );
  }

private:
  QToolButton		*but[15];
  CGAL::Qt_widget	*widget;
  QMainWindow		*window;
  int			nr_of_buttons;
  Input_mode            input_mode;

  Voronoi_diagram_layer<SDG_2>             *showVD;
  Skeleton_layer<SDG_2>                    *showSK;
  Sites_layer<SDG_2>                       *showSI;

};//end class

#endif // QT_LAYERS_TOOLBAR_H
