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

#ifndef CGAL_QT_LAYER_SHOW_CONSTRAINEDS_H
#define CGAL_QT_LAYER_SHOW_CONSTRAINEDS_H

#include <CGAL/IO/Qt_widget_layer.h>

namespace CGAL {

template <class T>
class Qt_layer_show_constraineds : public Qt_widget_layer {
public:
  typedef typename T::Edge            Edge;
  typedef typename T::Finite_edges_iterator
                                      Finite_edges_iterator;

  Qt_layer_show_constraineds(T &t) : tr(t){};

  void draw()
  {  
    Finite_edges_iterator it = tr.finite_edges_begin();
    *widget << CGAL::RED << CGAL::LineWidth(2);
    while(it != tr.finite_edges_end()) {
      if(tr.is_constrained(*it))
        *widget << tr.segment(it);
      ++it;
    }
  };
private:
  T	&tr;
  
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_CONSTRAINEDS_H
