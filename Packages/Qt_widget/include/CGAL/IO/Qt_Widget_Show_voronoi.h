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
// file          : include/CGAL/IO/Qt_Window_Show_voronoi.h
// package       : QT_window
// author(s)     : Ursu Radu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WINDOW_SHOW_VORONOI_H
#define CGAL_QT_WINDOW_SHOW_VORONOI_H


#include <qobject.h>
#include <CGAL/IO/Qt_Window.h>

namespace CGAL {

template <class T>
class Qt_widget_show_voronoi : public QObject
{
    //Q_OBJECT
public:
	Qt_widget_show_voronoi(){};

	void draw_view(T &t, Qt_widget &widget)
	{
		widget << CGAL::RED ;
		t.draw_dual(widget);
	};
	

};//end class 

} // namespace CGAL

#endif // CGAL_QT_WINDOW_GET_SEGMENT_H
