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
// author(s)     : Ursu Radu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_Widget.h>
#include "Qt_Widget_toolbar.h"

// icons

#include <CGAL/IO/pixmaps/handtool.xpm>
#include <CGAL/IO/pixmaps/zoom_in_rect.xpm>
#include <CGAL/IO/pixmaps/movepoint.xpm>
#include <CGAL/IO/pixmaps/no_tool.xpm>
#include <CGAL/IO/pixmaps/zoom_out.xpm>
#include <CGAL/IO/pixmaps/zoom_in.xpm>
#include <CGAL/IO/pixmaps/zoom_reg.xpm>
#include <CGAL/IO/pixmaps/polygon.xpm>


namespace CGAL {
  Tools_toolbar::Tools_toolbar(Qt_widget *w, QMainWindow *mw)
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
    but[0] = new QToolButton(QPixmap( (const char**)zoomin_reg_xpm ),
			     "Zoom in a certain region", 
			     0, 
			     this, 
			     SLOT(toolregion()), 
			     maintoolbar, 
			     "Zoom in a certain region");
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
		
    but[3] = new QToolButton(QPixmap( (const char**)zoomin_rect_xpm ),
			     "View rectangle", 
			     0, 
			     this, 
			     SLOT(zoominrect()), 
			     maintoolbar, 
			     "View rectangle");
		
    but[4] = new QToolButton(QPixmap( (const char**)notool_xpm ),
			     "Detach current tool", 
			     0, 
			     this, 
			     SLOT(notool()), 
			     maintoolbar, 
			     "Detach current tool");
		
    but[5] = new QToolButton(QPixmap( (const char**)handtool_xpm ),
						"Hand tool", 
						0, 
						this, 
						SLOT(handtool()), 
						maintoolbar, 
						"Hand tool");

    but[7] = new QToolButton(QPixmap( (const char**)polygon_xpm ),
						"Get Simple Polygon", 
						0, 
						this, 
						SLOT(get_s_polygon()), 
						maintoolbar, 
						"Get Simple Polygon");

		
    but[0]->setToggleButton(TRUE);
    but[3]->setToggleButton(TRUE);
    but[5]->setToggleButton(TRUE);
    but[7]->setToggleButton(TRUE);

    nr_of_buttons = 8;
  };

	
  //the definition of the slots
	
  
  void Tools_toolbar::toolregion()
  {
    if (but[0]->isOn())
    {
      if(is_active)
	but[activebutton]->toggle();
      *widget >> zoombut;
      activebutton = 0;
      is_active = TRUE;
    }
    else
    {
      widget->detach_current_tool();
      is_active = FALSE;
    }
  }

  void Tools_toolbar::zoominrect()
  {
    if (but[3]->isOn())
    {
	if(is_active)
	  but[activebutton]->toggle();
	*widget >> zoomrectbut;
	activebutton = 3;
	is_active = TRUE;
    }
    else
    {
	widget->detach_current_tool();
	is_active = FALSE;
    }
  }

  void Tools_toolbar::zoomin()
  {
    widget->zoom_in(2);
    widget->redraw();
  }
  void Tools_toolbar::zoomout()
  {
    widget->zoom_out(2);
    widget->redraw();
  }
  void Tools_toolbar::notool()
  {
    if(is_active)
      but[activebutton]->toggle();
    widget->detach_current_tool();
    is_active = FALSE;
  }
  

  void Tools_toolbar::handtool()
  {
    if (but[5]->isOn())
    {
      if(is_active)
	but[activebutton]->toggle();
      *widget >> handtoolbut;
      activebutton = 5;
      is_active = TRUE;
    } else {
	widget->detach_current_tool();
	is_active = FALSE;
    }
  }

  void Tools_toolbar::get_s_polygon()
  {
    if (but[7]->isOn())
    {
      if(is_active)
	but[activebutton]->toggle();
      *widget >> getsimplebut;
      activebutton = 7;
      is_active = TRUE;
    } else {
      widget->detach_current_tool();
      is_active = FALSE;
    }
  }


}//end namespace

#include "Qt_Widget_toolbar.moc"

#endif
