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
#include <CGAL/IO/Qt_layer_show_mouse_coordinates.h>
#include <CGAL/IO/Qt_layer_show_polygon.h>
#include <CGAL/IO/Qt_layer_show_greene_approximation.h>
#include <CGAL/IO/Qt_layer_show_ymonotone.h>
#include <CGAL/IO/Qt_layer_show_optimal_convex_partition.h>
#include <CGAL/IO/Qt_layer_show_polygon_points.h>

// icons
#include <CGAL/IO/pixmaps/mouse_coord.xpm>
#include <CGAL/IO/pixmaps/ymonotone.xpm>
#include <CGAL/IO/pixmaps/greene_approx.xpm>
#include <CGAL/IO/pixmaps/show_polygon.xpm>
#include <CGAL/IO/pixmaps/optimal_convex.xpm>
#include <CGAL/IO/pixmaps/points.xpm>


namespace CGAL {
  Layers_toolbar::Layers_toolbar(Qt_widget *w, QMainWindow *mw, Polygon *p) : 
     nr_of_buttons(0)
  {
    showMC = new Qt_layer_mouse_coordinates(*mw);
    showP = new Qt_layer_show_polygon<Polygon>(*p);
    showGA = new Qt_layer_show_greene_approx<Polygon>(*p);
    showYM = new Qt_layer_show_ymonotone<Polygon>(*p);
    showOC = new Qt_layer_show_optimal_convex<Polygon>(*p);
    showPP = new Qt_layer_show_polygon_points<Polygon>(*p);

    //set the widget
    widget = w;
    window = mw;
    window->statusBar();

    widget->attach(showMC);
    widget->attach(showP);
    widget->attach(showGA);
    widget->attach(showYM);
    widget->attach(showOC);
    widget->attach(showPP);

    maintoolbar = new QToolBar("tools", mw, QMainWindow::Top, TRUE, "Tools");
		
    but[0] = new QToolButton(QPixmap( (const char**)mouse_coord_xpm ),
					"Show Mouse Coordinates", 
					0, 
					this, 
					SLOT(show_coordinates()), 
					maintoolbar, 
					"Show Mouse Coordinates");
		
    
    but[1] = new QToolButton(QPixmap( (const char**)show_polygon_xpm ),
					"Show Polygon", 
					0, 
					this, 
					SLOT(show_polygon()), 
					maintoolbar, 
					"Show Polygon");

    
    but[2] = new QToolButton(QPixmap( (const char**)green_approx_xpm ),
					"Show Greene Approximation", 
					0, 
					this, 
					SLOT(show_greene_approx()), 
					maintoolbar, 
					"Show Greene Approximation");

    but[3] = new QToolButton(QPixmap( (const char**)ymonotone_xpm ),
					"Show Y Monotone Partition", 
					0, 
					this, 
					SLOT(show_ymonotone()), 
					maintoolbar, 
					"Show Y Monotone Partition");

    but[4] = new QToolButton(QPixmap( (const char**)optimal_convex_xpm ),
					"Show Optimal Convex Partition", 
					0, 
					this, 
					SLOT(show_optimal()), 
					maintoolbar, 
					"Show Optimal Convex Partition");

    but[5] = new QToolButton(QPixmap( (const char**)points_xpm ),
					"Show Polygon Vertices", 
					0, 
					this, 
					SLOT(show_points()), 
					maintoolbar, 
					"Show Polygon Vertices");


    
    nr_of_buttons = 6;		
    for(int i =0; i<nr_of_buttons; i++){
	but[i]->setToggleButton(TRUE);
	but[i]->toggle();
    }
			    
  }	
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
  
  void Layers_toolbar::show_polygon()
  {
    if (but[1]->isOn())
    {
      widget->activate(showP);
    } else {
      widget->deactivate(showP);
    }
    widget->redraw();
  }

  void Layers_toolbar::show_greene_approx()
  {
    if (but[2]->isOn())
    {
      widget->activate(showGA);
    } else {
      widget->deactivate(showGA);
    }
    widget->redraw();
  }
  void Layers_toolbar::show_ymonotone()
  {
    if (but[3]->isOn())
    {
      widget->activate(showYM);
    } else {
      widget->deactivate(showYM);
    }
    widget->redraw();
  }
  void Layers_toolbar::show_optimal()
  {
    if (but[4]->isOn())
    {
      widget->activate(showOC);
    } else {
      widget->deactivate(showOC);
    }
    widget->redraw();
  }
  void Layers_toolbar::show_points()
  {
    if (but[5]->isOn())
    {
      widget->activate(showPP);
    } else {
      widget->deactivate(showPP);
    }
    widget->redraw();
  }

  
}//end namespace

#include "Qt_widget_toolbar_layers.moc"

#endif
