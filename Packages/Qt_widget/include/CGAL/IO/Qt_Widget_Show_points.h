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
// file          : include/CGAL/IO/Qt_Window_Show_points.h
// package       : QT_window
// author(s)     : Ursu Radu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WINDOW_SHOW_POINTS_H
#define CGAL_QT_WINDOW_SHOW_POINTS_H

#include <CGAL/IO/Qt_Window.h>
#include <qobject.h>


namespace CGAL {

template <class T>
class Qt_widget_show_points : public QObject
{
    //Q_OBJECT
public:
	typedef typename T::Point						Point;
	typedef typename T::Segment					Segment;
	typedef typename T::Vertex					Vertex;
	typedef typename T::Vertex_iterator	Vertex_iterator;


	Qt_widget_show_points(){};


	void draw_view(T t, Qt_widget &widget)
	{
		Vertex v;
		Vertex_iterator it = t.vertices_begin(), 
						beyond = t.vertices_end();
		widget << CGAL::GREEN << CGAL::PointSize (3) << CGAL::PointStyle (CGAL::DISC);
		while(it != beyond)
		{
			v = *it;
			widget << v.point();
			++it;
		}
	};
	

};//end class 

} // namespace CGAL

#endif // CGAL_QT_WINDOW_GET_SEGMENT_H
