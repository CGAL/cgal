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
// file          : src/Qt_Window_toolbar_views.C
// package       : QT_window
// author(s)     : Ursu Radu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


#ifdef CGAL_USE_QT

#include "Qt_widget_toolbar_views.h"

// icons
#include <CGAL/IO/pixmaps/points.xpm>
#include <CGAL/IO/pixmaps/nearest_vertex.xpm>
#include <CGAL/IO/pixmaps/voronoi.xpm>
#include <CGAL/IO/pixmaps/triangulation.xpm>
#include <CGAL/IO/pixmaps/mouse_coord.xpm>


namespace CGAL {
  Views_toolbar::Views_toolbar(Qt_widget *w, QMainWindow *mw, Delaunay *t) : 
    dt(t), nr_of_buttons(0)
  {
    showT   = new Qt_view_show_triangulation< Delaunay >(*t);
    showV   = new Qt_view_show_voronoi< Delaunay >(*t);
    showP   = new Qt_view_show_points< Delaunay >(*t);
    showNV  = new Qt_view_nearest_vertex< Delaunay >(*t);
    showMC  = new Qt_view_mouse_coordinates(*mw);

    //set the widget
    widget = w;
    window = mw;
    window->statusBar();

    widget->attach(showT);
    widget->attach(showV);
    widget->attach(showNV);
    widget->attach(showP);
    widget->attach(showMC);
    widget->deactivate(showNV);
    

    maintoolbar = new QToolBar("tools", mw, QMainWindow::Top, TRUE, "Tools");
		

    but[0] = new QToolButton(QPixmap( (const char**)triangulation_xpm ),
			     "Show triangulation", 
			     0, 
			     this, 
			     SLOT(draw_triangulation()), 
			     maintoolbar, 
			     "Show triangulation");
		
    but[1] = new QToolButton(QPixmap( (const char**)voronoi_xpm ),
			     "Show Voronoi Diagram", 
			     0, 
			     this, 
			     SLOT(draw_voronoi()), 
			     maintoolbar, 
			     "Show Voronoi Diagram");
		
    but[2] = new QToolButton(QPixmap( (const char**)nearest_vertex_xpm ),
			     "Show Nearest Vertex", 
			     0, 
			     this, 
			     SLOT(draw_nearest_vertex()), 
			     maintoolbar, 
			     "Show Nearest Vertex");

    but[3] = new QToolButton(QPixmap( (const char**)points_xpm ),
				"Show Triangulation Points", 
				0, 
				this, 
				SLOT(draw_points()), 
				maintoolbar, 
				"Show Triangulation Points");
		
    but[4] = new QToolButton(QPixmap( (const char**)mouse_coord_xpm ),
				"Show Mouse Coordinates", 
				0, 
				this, 
				SLOT(show_coordinates()), 
				maintoolbar, 
				"Show Mouse Coordinates");
		
    nr_of_buttons = 5;
	
    for(int i =0; i<nr_of_buttons; i++)
    {
	but[i]->setToggleButton(TRUE);
	but[i]->toggle();
    }
    but[2]->toggle();
  }
  void Views_toolbar::draw_triangulation()
  {
    if (but[0]->isOn())
    {
      widget->activate(showT);
    } else {
      widget->deactivate(showT);
    }
    widget->redraw();
  }

  void Views_toolbar::draw_voronoi()
  {
    if (but[1]->isOn())
    {
      widget->activate(showV);
    } else {
      widget->deactivate(showV);
    }
    widget->redraw();
  }
  void Views_toolbar::draw_nearest_vertex()
  {
    if (but[2]->isOn())
    {
      widget->activate(showNV);
    } else {
      widget->deactivate(showNV);
    }
  }
  void Views_toolbar::draw_points()
  {
    if (but[3]->isOn())
    {
      widget->activate(showP);
    } else {
      widget->deactivate(showP);
    }
    widget->redraw();
  }
	
  void Views_toolbar::show_coordinates()
  {
    if (but[4]->isOn())
    {
      widget->activate(showMC);
      window->statusBar();
    } else {
      widget->deactivate(showMC);
      window->statusBar()->clear();
    }
  }

  
}//end namespace

#include "Qt_widget_toolbar_views.moc"

#endif
