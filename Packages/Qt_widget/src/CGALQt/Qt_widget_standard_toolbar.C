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
// file          : src/Qt_widget_standard_toolbar.C
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
#include <CGAL/IO/Qt_widget_standard_toolbar.h>

// icons
#include <CGAL/IO/pixmaps/zoom_in_rect.xpm>
#include <CGAL/IO/pixmaps/zoom_out.xpm>
#include <CGAL/IO/pixmaps/zoom_in.xpm>
#include <CGAL/IO/pixmaps/focus.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>



namespace CGAL {
  Qt_widget_standard_toolbar::Qt_widget_standard_toolbar(
	  Qt_widget *w, QMainWindow *mw)
  {

    w->attach_standard(&zoombut);
    w->attach_standard(&zoomrectbut);
    w->attach_standard(&handtoolbut);
    zoombut.deactivate();
    zoomrectbut.deactivate();
    handtoolbut.deactivate();
    //set the widget
    widget = w;

#if QT_VERSION < 300
  // for Qt 2.3 and before
  maintoolbar = new QToolBar("tools", mw, QMainWindow::Top, TRUE, "Tools");
#else
  // from Qt 3.0
  maintoolbar = new QToolBar(mw, "Tools");
  mw->addDockWindow (maintoolbar, "tools", DockTop, TRUE );
#endif
		
  but[0] = new QToolButton(maintoolbar, "notool");
  but[0]->setPixmap(QPixmap( (const char**)arrow_xpm ));
  
  but[1] = new QToolButton(QPixmap( (const char**)zoomin_xpm ),
			     "Zoom in", 
			     0, 
			     this, 
			     SLOT(zoomin()), 
			     maintoolbar, 
			     "Zoom in");
  but[2] = new QToolButton(QPixmap( (const char**)zoomout_xpm ),
			     "Zoom out", 
			     0, 
			     this, 
			     SLOT(zoomout()), 
			     maintoolbar, 
			     "Zoom out");

  but[3] = new QToolButton(maintoolbar, "focus");
  but[3]->setPixmap(QPixmap( (const char**)focus_xpm ));
  but[4] = new QToolButton(maintoolbar, "focus on region");
  but[4]->setPixmap(QPixmap( (const char**)zoomin_rect_xpm ));
  but[5] = new QToolButton(maintoolbar, "handtool");
  but[5]->setPixmap(QPixmap( (const char**)hand_xpm ));
    

	
    button_group = new QButtonGroup(0, "My_group");
    nr_of_buttons = 6;
    for(int i = 3; i<nr_of_buttons; i++){
      but[i]->setToggleButton(true);
      button_group->insert(but[i]);
    }
    but[0]->setToggleButton(true);
    button_group->insert(but[0]);
    button_group->setExclusive(true);
  
    connect(but[3], SIGNAL(stateChanged(int)),
        &zoombut, SLOT(stateChanged(int)));
    connect(but[4], SIGNAL(stateChanged(int)),
        &zoomrectbut, SLOT(stateChanged(int)));	
    connect(but[5], SIGNAL(stateChanged(int)),
        &handtoolbut, SLOT(stateChanged(int)));
  };
    void Qt_widget_standard_toolbar::zoomin()
  {
    widget->zoom_in(2);
    widget->redraw();
  }
  void Qt_widget_standard_toolbar::zoomout()
  {
    widget->zoom_out(2);
    widget->redraw();
  }


}//end namespace
#include "Qt_widget_standard_toolbar.moc"

#endif
