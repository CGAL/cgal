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
// file          : include/CGAL/IO/Qt_layer_show_triangulation.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_LAYER_SHOW_TRIANGULATION_H
#define CGAL_QT_LAYER_SHOW_TRIANGULATION_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_Triangulation_2.h>

namespace CGAL {

template <class T>
class Qt_layer_show_triangulation : public Qt_widget_layer
{
public:
	
  Qt_layer_show_triangulation(T &t) : tr(t){};


  void draw()
  {
    *widget << CGAL::BLUE; 
    *widget << tr;
  };
	
private:
  T &tr;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_TRIANGULATION_H
