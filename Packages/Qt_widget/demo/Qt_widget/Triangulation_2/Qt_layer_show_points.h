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
// file          : include/CGAL/IO/Qt_layer_show_points.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_LAYER_SHOW_POINTS_H
#define CGAL_QT_LAYER_SHOW_POINTS_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <qobject.h>


namespace CGAL {

template <class T>
class Qt_layer_show_points : public Qt_widget_layer {
public:
  typedef typename T::Point           Point;
  typedef typename T::Segment         Segment;
  typedef typename T::Vertex          Vertex;
  typedef typename T::Vertex_iterator	Vertex_iterator;

  Qt_layer_show_points(T &t) : tr(t){};

  void draw()
  {  
    Vertex_iterator it = tr.vertices_begin(), 
		beyond = tr.vertices_end();
    *widget << CGAL::GREEN << CGAL::PointSize (3) 
		<< CGAL::PointStyle (CGAL::DISC);    
    while(it != beyond) {      
      *widget << (*it).point();
      ++it;
    }
  };
private:
  T	&tr;
  
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_POINTS_H
