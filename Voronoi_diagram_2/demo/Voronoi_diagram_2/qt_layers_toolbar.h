// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

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
#include <CGAL/IO/pixmaps/circle.xpm>
#include <CGAL/IO/pixmaps/voronoi.xpm>
#include <CGAL/IO/pixmaps/mouse_coord.xpm>
#include <CGAL/IO/pixmaps/notool.xpm>
#include <CGAL/IO/pixmaps/holddown.xpm>
//#include "remove.xpm"

typedef enum { VD_POINT, VD_CIRCLE } Input_mode;

class Layers_toolbar : public QToolBar
{
  Q_OBJECT
public:
  Layers_toolbar(CGAL::Qt_widget *widget, VVD2* vvd,
		 const QString& label, QMainWindow* mainWindow,
		 QWidget* parent, bool newLine = FALSE,
		 const char* name = 0, WFlags f = 0 )
    : QToolBar(label, mainWindow, parent, newLine, name, f),
      nr_of_buttons(0), input_mode(VD_POINT)
  {
    showVD = new Voronoi_diagram_layer<VVD2>(vvd);
    showSI = new Sites_layer<VVD2>(vvd);

    // set the widget
    this->widget = widget;
    window = mainWindow;
    window->statusBar();

    widget->attach(showVD);
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

    but[2] = new QToolButton(QPixmap( (const char**)point_xpm ),
			     "Input points",
			     0,
			     this,
			     SLOT(input_points_mode()),
			     this,
			     "Input points");


    but[3] = new QToolButton(QPixmap( (const char**)circle_xpm ),
			     "Input circles",
			     0,
			     this,
			     SLOT(input_circles_mode()),
			     this,
			     "Input circles");
#if 1
    but[4] = new QToolButton(QPixmap( (const char**)notool_xpm ),
			     "Remove site",
			     0,
			     this,
			     SLOT(remove_mode()),
			     this,
			     "Remove site");
#else
    but[4] = new QToolButton(QPixmap( (const char**)remove_xpm ),
			     "Remove site",
			     0,
			     this,
			     SLOT(remove_mode()),
			     this,
			     "Remove site");
#endif
    but[5] = new QToolButton(QPixmap( (const char**)holddown_xpm ),
			     "Locate mode",
			     0,
			     this,
			     SLOT(locate_mode()),
			     this,
			     "Locate mode");

    nr_of_buttons = 6;
    for(int i = 0; i < nr_of_buttons; i++){
      but[i]->setToggleButton(TRUE);
    }

    but[0]->toggle();
    but[1]->toggle();
    but[2]->toggle();
  }

  ~Layers_toolbar() {
    delete showVD;
    delete showSI;
  }

  inline QToolBar* toolbar() { return this; };

  void set_vvd(VVD2* vvd) {
    showVD->set(vvd);
    showSI->set(vvd);
  }

signals:
  void new_object(CGAL::Object);
  void inputModeChanged(Input_mode);
  void insertModeChanged(bool);
  void locateModeChanged(bool);
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
    }

    widget->redraw();
  }

  void input_points_mode() {
    if ( !but[2]->isOn() ) {
      but[2]->toggle();
      return;
    }
    if ( input_mode == VD_CIRCLE ) { but[3]->toggle(); }

    input_mode = VD_POINT;
    emit inputModeChanged( VD_POINT );
  }

  void input_circles_mode() {
    if ( !but[3]->isOn() ) {
      but[3]->toggle();
      return;
    }
    if ( input_mode == VD_POINT ) { but[2]->toggle(); }

    input_mode = VD_CIRCLE;
    emit inputModeChanged( VD_CIRCLE );
  }

  void remove_mode() {
    emit insertModeChanged( but[4]->isOn() );
  }

  void locate_mode() {
    emit locateModeChanged( but[5]->isOn() );
  }

private:
  QToolButton		*but[15];
  CGAL::Qt_widget	*widget;
  QMainWindow		*window;
  int			nr_of_buttons;
  Input_mode            input_mode;

  Voronoi_diagram_layer<VVD2>             *showVD;
  Sites_layer<VVD2>                       *showSI;

};//end class

#endif // QT_LAYERS_TOOLBAR_H
