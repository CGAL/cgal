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

#ifndef CGAL_QT_WINDOW_NEAREST_VERTEX_H
#define CGAL_QT_WINDOW_NEAREST_VERTEX_H

#include <CGAL/IO/Qt_Window.h>
#include <qobject.h>




namespace CGAL {

template <class R, class T>
class Qt_widget_nearest_vertex : public QObject
{
    //Q_OBJECT
public:
	typedef typename T::Point							Point;
	typedef typename T::Segment						Segment;
	typedef typename T::Face_handle				Face_handle;
	typedef typename T::Vertex_handle			Vertex_handle;
	typedef typename T::Geom_traits::FT		FT;

	Point													oldPoint;
	bool													first_time;

	Qt_widget_nearest_vertex() : first_time(TRUE){};

	void widget_repainted(){
		first_time = TRUE;
	};

	void remove_view(Qt_widget &widget)
	{
		widget.lock();
		RasterOp old = widget.rasterOp();	//save the initial raster mode
		widget.setRasterOp(XorROP);
		widget << CGAL::GREEN << CGAL::PointSize (10)<< CGAL::PointStyle (CGAL::CIRCLE);
		widget << oldPoint;	
		widget.unlock();
		widget.setRasterOp(old);
		first_time = TRUE;
	}

	void mouseMove(QMouseEvent* e, T &t, Qt_widget &widget)
	{
		if (t.dimension()<1) return;
		FT
			x=static_cast<FT>(widget.x_real(e->x())),
			y=static_cast<FT>(widget.y_real(e->y()));

			Point p(x, y);
			RasterOp old = widget.rasterOp();	//save the initial raster mode
			widget.setRasterOp(XorROP);
			widget.lock();
			widget << Point(10, 10);
			Vertex_handle v = t.nearest_vertex(p);
			widget << CGAL::GREEN << CGAL::PointSize (10)<< CGAL::PointStyle (CGAL::CIRCLE);
			if(!first_time)
				widget << oldPoint;	
			widget << v->point();
			widget.unlock();
			widget.setRasterOp(old);
			oldPoint = v->point();
			first_time = FALSE;
	};
	
};//end class 

} // namespace CGAL

#endif // CGAL_QT_WINDOW_GET_SEGMENT_H
