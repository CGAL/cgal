// ======================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// file          : src/CGALQt/Qt_widget_standard_toolbar.C
// package       : Qt_widget (1.2.46)
// maintainer    : Laurent Rineau <rineau@clipper.ens.fr>
// author(s)     : Radu Ursu
// release       : $CGAL_Revision: CGAL-2.5-I-35 $
// release_date  : $CGAL_Date: 2002/10/18 $
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ======================================================================

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>

// icons
#include <CGAL/IO/pixmaps/zoom_in_rect.xpm>
#include <CGAL/IO/pixmaps/zoom_out.xpm>
#include <CGAL/IO/pixmaps/zoom_in.xpm>
#include <CGAL/IO/pixmaps/focus.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>
#include <CGAL/IO/pixmaps/back.xpm>
#include <CGAL/IO/pixmaps/forward.xpm>
#include <CGAL/IO/pixmaps/mouse_coord.xpm>

#include <qbuttongroup.h>
#include <qiconset.h>
#include <qmainwindow.h>
#include <qtoolbutton.h>
#include <qstring.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_focus.h>
#include <CGAL/IO/Qt_widget_zoomrect.h>
#include <CGAL/IO/Qt_widget_handtool.h>
#include <CGAL/IO/Qt_widget_show_mouse_coordinates.h>

namespace CGAL {
  Qt_widget_standard_toolbar::
  Qt_widget_standard_toolbar(Qt_widget *w, QMainWindow *mw,
#if QT_VERSION < 300
			     // for Qt 2.3 and before
			     QMainWindow::ToolBarDock dock,
#else
			     // from Qt 3.0
			     Dock dock,
#endif
			     bool newLine,
			     const char* name) :
    QToolBar("Qt_widget standard toolbar", mw)
  {
    setName(name);

    Qt_widget_focus* focusbut = 
      new Qt_widget_focus();
    Qt_widget_zoomrect* zoomrectbut = 
      new Qt_widget_zoomrect();
    Qt_widget_handtool* handtoolbut = 
      new Qt_widget_handtool();
    Qt_widget_show_mouse_coordinates* show_coord = 
      new Qt_widget_show_mouse_coordinates(*mw);
    // FIXME: thse objects should be destroyed in a destructor

    w->attach_standard(focusbut);
    w->attach_standard(zoomrectbut);
    w->attach_standard(handtoolbut);
    w->attach_standard(show_coord);
    show_coord->does_eat_events = false;
    focusbut->deactivate();
    zoomrectbut->deactivate();
    handtoolbut->deactivate();
    //set the widget
    widget = w;
    mw->statusBar();

#if QT_VERSION < 300
    // for Qt 2.3 and before
    mw->addToolBar(this, dock, newLine);
#else
    // from Qt 3.0
    mw->addDockWindow (this, dock, newLine);
#endif

    QIconSet set0(QPixmap( (const char**)arrow_small_xpm ),
                  QPixmap( (const char**)arrow_xpm ));
    QIconSet set1(QPixmap( (const char**)back_small_xpm ),
                  QPixmap( (const char**)back_xpm ));
    QIconSet set2(QPixmap( (const char**)forward_small_xpm ),
                  QPixmap( (const char**)forward_xpm ));

    QToolButton     *but[10];
    but[0] = new QToolButton(this, "nolayer");
    but[0]->setIconSet(set0);
    but[0]->setTextLabel("Deactivate Standard Layer");
  
    addSeparator();

    but[1] = new QToolButton(this, "History Back");
    but[1]->setIconSet(set1);
    but[1]->setTextLabel("History Back");

    but[2] = new QToolButton(this, "History Forward");
    but[2]->setIconSet(set2);
    but[2]->setTextLabel("History Forward");
  
    addSeparator();

    but[3] = new QToolButton(QPixmap( (const char**)zoomin_xpm ),
			     "Zoom in", 
			     0, 
			     this, 
			     SLOT(zoomin()), 
			     this, 
			     "Zoom in");
    but[3]->setTextLabel("Scaling factor X2");
    but[4] = new QToolButton(QPixmap( (const char**)zoomout_xpm ),
			     "Zoom out", 
			     0, 
			     this, 
			     SLOT(zoomout()), 
			     this, 
			     "Zoom out");
    but[4]->setTextLabel("Scaling factor 1/2");
  
    but[5] = new QToolButton(this, "focus on region");
    but[5]->setPixmap(QPixmap( (const char**)zoomin_rect_xpm ));
    but[5]->setTextLabel("Focus on region");

    addSeparator();

    but[6] = new QToolButton(this, "focus");
    but[6]->setPixmap(QPixmap( (const char**)focus_xpm ));
    but[6]->setTextLabel("Focus on point");
    but[7] = new QToolButton(this, "handtool");
    but[7]->setPixmap(QPixmap( (const char**)hand_xpm ));
    but[7]->setTextLabel("Pan tool");

    addSeparator();

    but[8] = new QToolButton(this, "mouse");
    but[8]->setPixmap(QPixmap( (const char**)mouse_coord_xpm) );
    but[8]->setTextLabel("Mouse Coordinates");

    QButtonGroup* button_group = new QButtonGroup(0, "My_group");
    // FIXME: this QButtonGroup has no parent and should be destroyed 
    int nr_of_buttons = 9;
    for(int i = 5; i<nr_of_buttons-1; i++){
      but[i]->setToggleButton(true);
      button_group->insert(but[i]);
    }
    but[0]->setToggleButton(true);
    but[8]->setToggleButton(true);
    but[8]->toggle();

    button_group->insert(but[0]);
    button_group->setExclusive(true);
  
    connect(but[1], SIGNAL(clicked()),
	    this, SLOT(back()));
    connect(but[2], SIGNAL(clicked()),
	    this, SLOT(forward()));
    connect(but[5], SIGNAL(stateChanged(int)),
        zoomrectbut, SLOT(stateChanged(int)));
    connect(but[6], SIGNAL(stateChanged(int)),
        focusbut, SLOT(stateChanged(int)));	
    connect(but[7], SIGNAL(stateChanged(int)),
        handtoolbut, SLOT(stateChanged(int)));
    connect(but[8], SIGNAL(stateChanged(int)),
        show_coord, SLOT(stateChanged(int)));
    
    connect(widget, SIGNAL(set_back_enabled(bool)),
            but[1], SLOT(setEnabled(bool)));
    connect(widget, SIGNAL(set_forward_enabled(bool)),
            but[2], SLOT(setEnabled(bool)));
    widget->clear_history();
  };
  void Qt_widget_standard_toolbar::zoomin()
  {
    widget->zoom(2); //
  }
  void Qt_widget_standard_toolbar::zoomout()
  {
    widget->zoom(0.5);
  }
  void Qt_widget_standard_toolbar::back()
  {
    widget->back();
  }
    void Qt_widget_standard_toolbar::forward()
  {
    widget->forth();
  }

}//end namespace
#include "Qt_widget_standard_toolbar.moc"

#endif
