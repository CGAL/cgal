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
#include <CGAL/IO/pixmaps/zoom_reg.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>



namespace CGAL {
  Qt_widget_standard_toolbar::Qt_widget_standard_toolbar(
	  Qt_widget *w, QMainWindow *mw)
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
		
    but[3] = new QToolButton(QPixmap( (const char**)zoomin_reg_xpm ),
			     "Zoom in a certain region", 
			     0, 
			     this, 
			     SLOT(toolregion()), 
			     maintoolbar, 
			     "Zoom in a certain region");
    
    but[4] = new QToolButton(QPixmap( (const char**)zoomin_rect_xpm ),
			     "Layer rectangle", 
			     0, 
			     this, 
			     SLOT(zoominrect()), 
			     maintoolbar, 
			     "Layer rectangle");
		
    		
    but[5] = new QToolButton(QPixmap( (const char**)hand_xpm ),
			"Hand tool", 
			0, 
			this, 
			SLOT(handtool()), 
			maintoolbar, 
			"Hand tool");		
		
    but[3]->setToggleButton(TRUE);
    but[4]->setToggleButton(TRUE);
    but[5]->setToggleButton(TRUE);
    nr_of_buttons = 6;

    connect(w, SIGNAL(detached_standard_tool()), this, SLOT(toggle_button()));
  };

  //the definition of the slots
  void Qt_widget_standard_toolbar::toggle_button()
  {
    if(is_active) {
      but[activebutton]->toggle();
      is_active = false;
    }
  }	

  void Qt_widget_standard_toolbar::toolregion()
  {
    if (but[3]->isOn())
    {
      widget->attach_standard(&zoombut);
      activebutton = 3;
      is_active = true;
    }
    else
    {
      is_active = false;
      widget->detach_current_standard_tool();
    }
  }

  void Qt_widget_standard_toolbar::zoominrect()
  {
    if (but[4]->isOn())
    {
      widget->attach_standard(&zoomrectbut);
      activebutton = 4;
      is_active = true;
    }
    else
    {
      is_active = false;
      widget->detach_current_standard_tool();
    }
  }

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
  void Qt_widget_standard_toolbar::notool()
  {
    if(is_active) {
      but[activebutton]->toggle();
      widget->detach_current_standard_tool();
      is_active = false;
    }
  }


  void Qt_widget_standard_toolbar::handtool()
  {
    if (but[5]->isOn())
    {      
      widget->attach_standard(&handtoolbut);
      activebutton = 5;
      is_active = true;
    } else {
      is_active = false;
      widget->detach_current_standard_tool();
    }
  }

}//end namespace
#include "Qt_widget_standard_toolbar.moc"

#endif
