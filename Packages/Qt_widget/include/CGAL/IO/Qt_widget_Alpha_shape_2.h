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
// file          : include/CGAL/IO/Qt_widget_Alpha_shape_2.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


#ifndef CGAL_QT_WIDGET_ALPHA_SHAPE_2_H
#define CGAL_QT_WIDGET_ALPHA_SHAPE_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Alpha_shape_2.h>
namespace CGAL{

template< class Dt >
Qt_widget&
operator << ( Qt_widget& ws, const CGAL::Alpha_shape_2<Dt>& As)
{
  //return As.op_window(ws);
  typedef typename Alpha_shape_2<Dt>::Alpha_shape_edges_iterator 
                    Edges_iterator;
  typedef typename Alpha_shape_2<Dt>::Segment Segment_2;
  if (As.get_mode() == Alpha_shape_2<Dt>::REGULARIZED) 
  { 
    for (Edges_iterator edge_alpha_it = As.alpha_shape_edges_begin();
         edge_alpha_it != As.alpha_shape_edges_end(); edge_alpha_it++)
    {
      ws << As.segment(*edge_alpha_it);
    }//endfor

  } else {
    for (Edges_iterator edge_alpha_it = As.alpha_shape_edges_begin();
         edge_alpha_it != As.alpha_shape_edges_end(); edge_alpha_it++)
    {
      ws << As.segment(*edge_alpha_it);
    }//endfor
  }
      return ws;
}

}//end namespace CGAL

#endif
