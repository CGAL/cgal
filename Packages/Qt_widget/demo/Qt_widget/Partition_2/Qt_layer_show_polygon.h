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
// file          : include/CGAL/IO/Qt_layer_show_polygon.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_LAYER_SHOW_POLYGON_H
#define CGAL_QT_LAYER_SHOW_POLYGON_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <qobject.h>




namespace CGAL {

template <class T>
class Qt_layer_show_polygon : public Qt_widget_layer
{
public:
  
  Qt_layer_show_polygon(T &p) : polygon(p){};
  void draw()
  {
    *widget << LineWidth(3);
    *widget << CGAL::BLUE; 
    *widget << polygon;
    *widget << LineWidth(1);
    *widget << CGAL::WHITE; 
    *widget << polygon;
  };
	
private:
  T &polygon;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_POLYGON_H
