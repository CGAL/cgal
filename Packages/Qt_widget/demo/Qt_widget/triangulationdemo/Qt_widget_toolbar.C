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
#include "Qt_widget_toolbar.h"

// icons
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/handtool.xpm>
#include <CGAL/IO/pixmaps/zoom_in_rect.xpm>
#include <CGAL/IO/pixmaps/movepoint.xpm>
#include <CGAL/IO/pixmaps/no_tool.xpm>
#include <CGAL/IO/pixmaps/zoom_out.xpm>
#include <CGAL/IO/pixmaps/zoom_in.xpm>
#include <CGAL/IO/pixmaps/zoom_reg.xpm>
#include <CGAL/IO/pixmaps/line.xpm>


namespace CGAL {
  Tools_toolbar::Tools_toolbar(Qt_widget *w, QMainWindow *mw, Delaunay *t) : dt(t)
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
		

    but[0] = new QToolButton(QPixmap( (const char**)point_xpm ),
			     "Point Tool", 
			     0, 
			     this, 
			     SLOT(pointtool()), 
			     maintoolbar, 
			     "Point Tool");

    but[1] = new QToolButton(QPixmap( (const char**)line_xpm ),
			     "Line Tool", 
			     0, 
			     this, 
			     SLOT(linetool()), 
			     maintoolbar, 
			     "Line Tool");
    but[2] = new QToolButton(QPixmap( (const char**)zoomin_reg_xpm ),
			     "Zoom in a certain region", 
			     0, 
			     this, 
			     SLOT(toolregion()), 
			     maintoolbar, 
			     "Zoom in a certain region");
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
		
    but[5] = new QToolButton(QPixmap( (const char**)zoomin_rect_xpm ),
			     "View rectangle", 
			     0, 
			     this, 
			     SLOT(zoominrect()), 
			     maintoolbar, 
			     "View rectangle");
		
    but[6] = new QToolButton(QPixmap( (const char**)notool_xpm ),
			     "Detach current tool", 
			     0, 
			     this, 
			     SLOT(notool()), 
			     maintoolbar, 
			     "Detach current tool");
		
    but[7] = new QToolButton(QPixmap( (const char**)handtool_xpm ),
			"Hand tool", 
			0, 
			this, 
			SLOT(handtool()), 
			maintoolbar, 
			"Hand tool");

    but[8] = new QToolButton(QPixmap( (const char**)movepoint_xpm ),
			     "Move current point", 
			     0, 
			     this, 
			     SLOT(movepoint()), 
			     maintoolbar, 
			     "Move point");
		
		
    but[0]->setToggleButton(TRUE);
    but[1]->setToggleButton(TRUE);
    but[2]->setToggleButton(TRUE);
    but[5]->setToggleButton(TRUE);
    but[7]->setToggleButton(TRUE);
    but[8]->setToggleButton(TRUE);

		nr_of_buttons = 9;
  };

      
	
  //the definition of the slots
	
  void Tools_toolbar::linetool()
  {
    if (but[1]->isOn())
    {
	
	if(is_active)
	  but[activebutton]->toggle();
	*widget >> linebut;
	activebutton = 1;
	is_active = TRUE;
    }
    else
    {
	widget->detach_current_tool();
	is_active = FALSE;
    }
  }
  void Tools_toolbar::pointtool()
  {
    if (but[0]->isOn())
    {
	
	if(is_active)
	  but[activebutton]->toggle();
	*widget >> pointbut;
	activebutton = 0;
	is_active = TRUE;
    }
    else
    {
	widget->detach_current_tool();
	is_active = FALSE;
    }
  }
  void Tools_toolbar::toolregion()
  {
    if (but[2]->isOn())
    {
	if(is_active)
	  but[activebutton]->toggle();
	*widget >> zoombut;
	activebutton = 2;
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
    if (but[5]->isOn())
    {	
	if(is_active)
	  but[activebutton]->toggle();
	*widget >> zoomrectbut;
	activebutton = 5;
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
  void Tools_toolbar::movepoint()
  {
    if (but[8]->isOn())
    {
			widget->detach_current_tool();
			if(is_active)
				but[activebutton]->toggle();
			movepointbut.set_Delaunay(dt);
			*widget >> movepointbut;
			activebutton = 8;
			is_active = TRUE;
    }
    else
    {
			widget->detach_current_tool();
			is_active = FALSE;
    }
  }


	void Tools_toolbar::handtool()
	{
		if (but[7]->isOn())
		{
			widget->detach_current_tool();
			if(is_active)
				but[activebutton]->toggle();
			*widget >> handtoolbut;
			activebutton = 7;
			is_active = TRUE;
		} else {
			widget->detach_current_tool();
			is_active = FALSE;
		}
	}

}//end namespace

#include "Qt_widget_toolbar.moc"

#endif
