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
// file          : include/CGAL/IO/Qt_Window_toolbar_views.h
// package       : QT_window
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WINDOW_TOOLBAR_VIEWS_H
#define CGAL_QT_WINDOW_TOOLBAR_VIEWS_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Qt_widget.h>

#include "min_quadrilateral_views.h"
#include <CGAL/IO/Qt_view_show_mouse_coordinates.h>

#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qstatusbar.h>

typedef double Coord_type;
typedef CGAL::Cartesian<Coord_type>  Rp;
typedef Rp::Point_2   Point;

namespace CGAL {

class Views_toolbar : public QObject
{
	Q_OBJECT
public:
  Views_toolbar(Qt_widget *w, QMainWindow *mw, std::list<Point>	*l_of_p);
	~Views_toolbar()
	{
		delete maintoolbar;
	};
	QToolBar*	toolbar(){return maintoolbar;};

signals:
	void new_object(CGAL::Object);
	
private slots:
	void show_coordinates();
	void show_points();
	void show_strip();
	void show_rectangle();
	void show_parallelogram();

private:
	QToolBar		*maintoolbar;
	QToolButton		*but[10];
	Qt_widget		*widget;
	QMainWindow		*window;
	int			nr_of_buttons;
	
	CGAL::Qt_view_mouse_coordinates	  *showMC;
	Qt_view_show_points<Rp>		  *showP;
	Qt_view_show_parallelogram<Rp>	  *showPL;
	Qt_view_show_strip<Rp>	  *showLS;
	Qt_view_show_rectangle<Rp>	  *showR;
};//end class

};//end namespace

#endif
