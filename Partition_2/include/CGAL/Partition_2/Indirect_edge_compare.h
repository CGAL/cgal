// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_INDIRECT_EDGE_COMPARE_H
#define CGAL_INDIRECT_EDGE_COMPARE_H

#include <CGAL/license/Partition_2.h>


namespace CGAL {

//
// given circulators to endpoints of two edges, sorts the edges that come
// next (in the direction of circulator) from right to left. This ordering 
// makes finding the edge directly left of a given edge (needed for the 
// y-monotone decomposition algorithm) easy. 
//
template <class ForwardCirculator, class Traits>
class Indirect_edge_compare 
{
   public:
     typedef typename Traits::Compare_y_2        Compare_y_2;
     typedef typename Traits::Compare_x_2        Compare_x_2;
     typedef typename Traits::Construct_line_2   Construct_line_2;
     typedef typename Traits::Compare_x_at_y_2   Compare_x_at_y_2;
     typedef typename Traits::Is_horizontal_2    Is_horizontal_2;
     typedef typename Traits::Line_2             Line_2;
     typedef typename Traits::Point_2            Point_2;

     Indirect_edge_compare() : 
          _compare_y_2(Traits().compare_y_2_object()),
          _compare_x_2(Traits().compare_x_2_object()),
          _construct_line_2(Traits().construct_line_2_object()),
          _compare_x_at_y_2(Traits().compare_x_at_y_2_object()),
          _is_horizontal_2(Traits().is_horizontal_2_object())
     { }
     
     // determines if the edge (edge_vtx_1, edge_vtx_1++) has a larger
     // x value than vertex.x at y-value vertex.y
     bool
     larger_x_at_vertex_y(ForwardCirculator edge_vtx_1, 
                          ForwardCirculator vertex) const
     {
        ForwardCirculator edge_vtx_2 = edge_vtx_1;
        edge_vtx_2++;
        // check for horizontal edge
        if (_compare_y_2((*edge_vtx_1), (*edge_vtx_2)) == EQUAL)  
        { 
            // compare the smaller x and vertex x
           if (_compare_x_2(*edge_vtx_1, *edge_vtx_2) == SMALLER)
              return _compare_x_2(*edge_vtx_1, *vertex) == LARGER;
           else
              return _compare_x_2(*edge_vtx_2, *vertex) == LARGER;
        }
        else 
        { 
           // construct supporting line for edge
           Line_2  line = _construct_line_2(*edge_vtx_1, *edge_vtx_2);
           return _compare_x_at_y_2(*vertex, line) == SMALLER;
        }
     }               

     bool 
     operator()(ForwardCirculator p, ForwardCirculator q) const
     {
        ForwardCirculator after_p = p;
        after_p++;
        ForwardCirculator after_q = q;
        after_q++;

        if (p == q && after_p == after_q) return false;

        if (p == after_q) 
          return larger_x_at_vertex_y(p, q);

        if (after_p == q) 
          return !larger_x_at_vertex_y(q, p);

        if (p == q) 
          return larger_x_at_vertex_y(p, after_q);

        if (after_p == after_q) 
          return larger_x_at_vertex_y(p, q);

        // else neither endpoint is shared
        // construct supporting line
        Line_2  l_p = _construct_line_2(*p, *after_p);
        if (_is_horizontal_2(l_p)) 
        {
            Line_2  l_q = _construct_line_2(*q, *after_q);

            if (_is_horizontal_2(l_q))  
            {                         
                 Point_2 p_max;
                 Point_2 q_max;
                 if (_compare_x_2(*p, *after_p) == SMALLER)
                    p_max = *after_p;
                 else
                    p_max = *p;
                 if (_compare_x_2(*q, *after_q) == SMALLER)
                    q_max = *after_q;
                 else
                    q_max = *q;
                 return (_compare_x_2(p_max, q_max) == LARGER);
            }
            else  // p and after_p must both be on same side of l_q
            {
                 return (_compare_x_at_y_2(*p, l_q) == LARGER);
            }
        }

        bool q_larger_x = _compare_x_at_y_2(*q, l_p) == SMALLER;
        bool after_q_larger_x = _compare_x_at_y_2(*after_q, l_p) == SMALLER;

        if (q_larger_x == after_q_larger_x)
            return q_larger_x;
        // else one smaller and one larger
        // construct the other line
        Line_2 l_q = _construct_line_2(*q, *after_q); 
        if (_is_horizontal_2(l_q))     // p is not horizontal
        {
            return _compare_x_at_y_2((*q), l_p) == LARGER;
        }
        return _compare_x_at_y_2((*p), l_q) != SMALLER;
     }

   private:
     Compare_y_2      _compare_y_2;
     Compare_x_2      _compare_x_2;
     Construct_line_2 _construct_line_2;
     Compare_x_at_y_2 _compare_x_at_y_2;
     Is_horizontal_2  _is_horizontal_2;
};

}



#endif // CGAL_INDIRECT_EDGE_COMPARE_H
