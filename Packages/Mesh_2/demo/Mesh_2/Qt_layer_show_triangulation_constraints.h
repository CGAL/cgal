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
// file          : include/CGAL/IO/Qt_layer_show_triangulation_constrints.h
// package       : Qt_widget
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_LAYER_SHOW_TRIANGULATION_CONSTRAINTS_H
#define CGAL_QT_LAYER_SHOW_TRIANGULATION_CONSTRAINTS_H

#include <CGAL/IO/Qt_widget_layer.h>

namespace CGAL {

template <class T>
class Qt_layer_show_triangulation_constraints : public Qt_widget_layer
{
public:
	
  Qt_layer_show_triangulation_constraints(T *&t) : tr(t){};

  void draw()
  {
    widget->lock();
    for(typename T::Edge_iterator it=tr->edges_begin();
	it!=tr->edges_end();
	it++)
      if(tr->is_constrained(*it))
	  *widget << CGAL::RED << tr->segment(*it);
    widget->unlock();  
  };
	
private:
  T *&tr;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_TRIANGULATION_CONSTRAINTS_H
