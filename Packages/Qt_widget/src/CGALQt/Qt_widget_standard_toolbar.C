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
#include <CGAL/IO/pixmaps/back.xpm>
#include <CGAL/IO/pixmaps/forward.xpm>



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
    widget->clear_history();

#if QT_VERSION < 300
  // for Qt 2.3 and before
  maintoolbar = new QToolBar("Qt_widget standard toolbar", mw, QMainWindow::Top, TRUE, "std_toolbar");
#else
  // from Qt 3.0
  maintoolbar = new QToolBar(mw, "std_toolbar");
  mw->addDockWindow (maintoolbar, "Qt_widget standard toolbar",
		     DockTop, TRUE );
#endif
		
  but[0] = new QToolButton(maintoolbar, "notool");
  but[0]->setPixmap(QPixmap( (const char**)arrow_xpm ));
  
  but[1] = new QToolButton(QPixmap( (const char**)back_xpm ),
			     "Back", 
			     0, 
			     this, 
			     SLOT(back()), 
			     maintoolbar, 
			     "Back");
  but[2] = new QToolButton(QPixmap( (const char**)forward_xpm ),
			     "Forward", 
			     0, 
			     this, 
			     SLOT(forward()), 
			     maintoolbar, 
			     "Forward");

  
  but[3] = new QToolButton(QPixmap( (const char**)zoomin_xpm ),
			     "Zoom in", 
			     0, 
			     this, 
			     SLOT(zoomin()), 
			     maintoolbar, 
			     "Zoom in");
  but[4] = new QToolButton(QPixmap( (const char**)zoomout_xpm ),
			     "Zoom out", 
			     0, 
			     this, 
			     SLOT(zoomout()), 
			     maintoolbar, 
			     "Zoom out");
  
  but[5] = new QToolButton(maintoolbar, "focus");
  but[5]->setPixmap(QPixmap( (const char**)focus_xpm ));
  but[6] = new QToolButton(maintoolbar, "focus on region");
  but[6]->setPixmap(QPixmap( (const char**)zoomin_rect_xpm ));
  but[7] = new QToolButton(maintoolbar, "handtool");
  but[7]->setPixmap(QPixmap( (const char**)hand_xpm ));
    

	
    button_group = new QButtonGroup(0, "My_group");
    nr_of_buttons = 8;
    for(int i = 5; i<nr_of_buttons; i++){
      but[i]->setToggleButton(true);
      button_group->insert(but[i]);
    }
    but[0]->setToggleButton(true);
    //but[1]->setEnabled(false);
    //but[2]->setEnabled(false);
    button_group->insert(but[0]);
    button_group->setExclusive(true);
  
    connect(but[5], SIGNAL(stateChanged(int)),
        &zoombut, SLOT(stateChanged(int)));
    connect(but[6], SIGNAL(stateChanged(int)),
        &zoomrectbut, SLOT(stateChanged(int)));	
    connect(but[7], SIGNAL(stateChanged(int)),
        &handtoolbut, SLOT(stateChanged(int)));
    
    connect(widget, SIGNAL(set_back_enabled(bool)),
            but[1], SLOT(setEnabled(bool)));
    connect(widget, SIGNAL(set_forward_enabled(bool)),
            but[2], SLOT(setEnabled(bool)));
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
  void Qt_widget_standard_toolbar::back()
  {
    widget->back();
    widget->redraw();
  }
    void Qt_widget_standard_toolbar::forward()
  {
    widget->forth();
    widget->redraw();
  }

}//end namespace
#include "Qt_widget_standard_toolbar.moc"

#endif
