// Copyright (c) 1997-2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Radu Ursu


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
