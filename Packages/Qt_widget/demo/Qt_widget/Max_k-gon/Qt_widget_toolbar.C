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
// file          : src/Qt_Window_toolbar.C
// package       : QT_window
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_widget.h>
#include "Qt_widget_toolbar.h"

// icons
#include <CGAL/IO/pixmaps/movepoint.xpm>
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>


namespace CGAL {
  Tools_toolbar::Tools_toolbar(Qt_widget *w, QMainWindow *mw, std::list<Point> &l1)
  {
    //when it is created, the toolbar has 0 buttons
    nr_of_buttons = 0;
    //set the widget
    widget = w;
    is_active = FALSE;

#if QT_VERSION < 300
  // for Qt 2.3 and before
  maintoolbar = new QToolBar("tools", mw, QMainWindow::Top, TRUE, "Tools");
#else
  // from Qt 3.0
  maintoolbar = new QToolBar(mw, "Tools");
  mw->addDockWindow (maintoolbar, "tools", DockTop, TRUE );
#endif
		

  but[0] = new QToolButton(QPixmap( (const char**)arrow_xpm ),
			     "Detach current tool", 
			     0, 
			     this, 
			     SLOT(notool()), 
			     maintoolbar, 
			     "Detach current tool");

  but[1] = new QToolButton(QPixmap( (const char**)point_xpm ),
			     "Point Tool", 
			     0, 
			     this, 
			     SLOT(pointtool()), 
			     maintoolbar, 
			     "Point Tool");

  
  but[1]->setToggleButton(TRUE);

  nr_of_buttons = 2;

  connect(w, SIGNAL(detached_tool()), this, SLOT(toggle_button()));
};

      
  //the definition of the slots
  void Tools_toolbar::toggle_button ()
  {
    if(is_active) {
      but[activebutton]->toggle();
      is_active = false;
    }
  }	
  
  void Tools_toolbar::pointtool()
  {
    if (but[1]->isOn())
    {
      widget->attach(pointbut);
      activebutton = 1;
      is_active = true;
    }
    else
    {
      is_active = false;
      widget->detach_current_tool();
    }
  }

  void Tools_toolbar::notool()
  {
    if(is_active) {
      widget->detach_current_tool();
      is_active = false;
    }
  }
  

}//end namespace

#include "Qt_widget_toolbar.moc"

#endif
