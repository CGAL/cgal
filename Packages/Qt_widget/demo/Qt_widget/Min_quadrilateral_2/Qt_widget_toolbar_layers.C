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
// file          : src/Qt_Window_toolbar_layers.C
// package       : QT_window
// author(s)     : Ursu Radu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


#ifdef CGAL_USE_QT

#include "Qt_widget_toolbar_layers.h"

// icons
#include <CGAL/IO/pixmaps/points.xpm>
#include <CGAL/IO/pixmaps/line.xpm>
#include <CGAL/IO/pixmaps/min_rectangle.xpm>
#include <CGAL/IO/pixmaps/min_parallelogram.xpm>
#include <CGAL/IO/pixmaps/mouse_coord.xpm>


namespace CGAL {
  Layers_toolbar::Layers_toolbar(Qt_widget *w, QMainWindow *mw, std::list<Point>	*l_of_p) : 
    widget(w), window(mw), nr_of_buttons(0)
  {
      
    showMC  = new Qt_layer_mouse_coordinates(*mw);
    showPL  = new Qt_layer_show_parallelogram<Rp>(l_of_p);
    showP   = new Qt_layer_show_points<Rp>(l_of_p);
    showLS  = new Qt_layer_show_strip<Rp>(l_of_p);
    showR   = new Qt_layer_show_rectangle<Rp>(l_of_p);
    
    //set the widget
    window->statusBar();

    widget->attach(showMC);
    widget->attach(showR);
    widget->attach(showPL);
    widget->attach(showP);
    widget->attach(showLS);

    maintoolbar = new QToolBar("tools", mw, QMainWindow::Top, TRUE, "Tools");
		
    but[0] = new QToolButton(QPixmap( (const char**)mouse_coord_xpm ),
				"Show Mouse Coordinates", 
				0, 
				this, 
				SLOT(show_coordinates()), 
				maintoolbar, 
				"Show Mouse Coordinates");
		
    but[1] = new QToolButton(QPixmap( (const char**)points_xpm ),
				"Show Points", 
				0, 
				this, 
				SLOT(show_points()), 
				maintoolbar, 
				"Show Points");

    but[2] = new QToolButton(QPixmap( (const char**)min_parallelogram_xpm ),
				"Show Parallelogram", 
				0, 
				this, 
				SLOT(show_parallelogram()), 
				maintoolbar, 
				"Show Parallelogram");

    but[3] = new QToolButton(QPixmap( (const char**)line_xpm ),
				"Show The Line Strip", 
				0, 
				this, 
				SLOT(show_strip()), 
				maintoolbar, 
				"Show The Line Strip");

    but[4] = new QToolButton(QPixmap( (const char**)min_rectangle_xpm ),
				"Show The Minimum Rectangle", 
				0, 
				this, 
				SLOT(show_rectangle()), 
				maintoolbar, 
				"Show The Minimum Rectangle");

    nr_of_buttons = 5;
	
    for(int i =0; i<nr_of_buttons; i++)
    {
	but[i]->setToggleButton(TRUE);
	but[i]->toggle();
    }
  }//end constructor
  	
  void Layers_toolbar::show_coordinates()
  {
    if (but[0]->isOn())
    {
      widget->activate(showMC);
      window->statusBar();
    } else {
      widget->deactivate(showMC);
      window->statusBar()->clear();
    }
  }
  void Layers_toolbar::show_points()
  {
    if (but[1]->isOn())
    {
      widget->activate(showP);
      widget->redraw();
    } else {
      widget->deactivate(showP);
      widget->redraw();
    }
  }

  void Layers_toolbar::show_parallelogram()
  {
    if (but[2]->isOn())
    {
      widget->activate(showPL);
      widget->redraw();
    } else {
      widget->deactivate(showPL);
      widget->redraw();
    }
  }
  void Layers_toolbar::show_strip()
  {
    if (but[3]->isOn())
    {
      widget->activate(showLS);
      widget->redraw();
    } else {
      widget->deactivate(showLS);
      widget->redraw();
    }
  }

  void Layers_toolbar::show_rectangle()
  {
    if (but[4]->isOn())
    {
      widget->activate(showR);
      widget->redraw();
    } else {
      widget->deactivate(showR);
      widget->redraw();
    }
  }


}//end namespace

#include "Qt_widget_toolbar_layers.moc"

#endif
