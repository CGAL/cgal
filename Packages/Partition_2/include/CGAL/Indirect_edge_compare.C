// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : include/CGAL/Indirect_edge_compare.C
// package       : $CGAL_Package: Partition_2 $
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
// chapter       : Planar Polygon Partitioning
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: Comparison of edges represented by circulators to endpoints.
// ============================================================================

#include <CGAL/Indirect_edge_compare.h>

namespace CGAL {

template <class Traits>
template <class ForwardCirculator>
bool
Indirect_edge_compare<Traits>:: larger_x_at_vertex_y(
                                             ForwardCirculator edge_vtx_1, 
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
      // compute x value at vertex's y value
// ??? CHANGE THIS when the function Compare_x_at_y is available
//         return _compare_x_at_y(*vertex, line) == SMALLER;
// ???
      return line.x_at_y((*vertex).y()) > (*vertex).x();
   }
}

template <class Traits>
template <class ForwardCirculator>
bool 
Indirect_edge_compare<Traits>::operator()(ForwardCirculator p, 
                                          ForwardCirculator q) const
{ 
   ForwardCirculator after_p = p;
   after_p++;
   ForwardCirculator after_q = q;
   after_q++;
   if (p == after_q) 
   {
     return larger_x_at_vertex_y(p, q);
   }
   else if (after_p == q) 
   {
     return !larger_x_at_vertex_y(q, p);
   }
   else if (p == q) 
   {
     return larger_x_at_vertex_y(p, after_q);
   }
   else if (after_p == after_q) 
   {
     return larger_x_at_vertex_y(p, q);
   }
   else // neither endpoint is shared
   {
     // construct supporting line
     Line_2  l_p = _construct_line_2(*p, *after_p);
     if (l_p.is_horizontal()) 
     {
         Line_2  l_q = _construct_line_2(*q, *after_q);
         if (l_q.is_horizontal())  // shouldn't ever happen, since these
         {                         // can't both be in sweep structure at
                                   // the same time
            return std::max((*p).x(), (*after_p).x()) > 
                   std::max((*q).x(), (*after_q).x());
         }
         else  // p and after_p must both be on same side of l_q
         {
            return (*p).x() > l_q.x_at_y((*p).y());
         }
     }
     else  
     {
        bool q_larger_x = l_p.x_at_y((*q).y()) > (*q).x();
        bool after_q_larger_x = l_p.x_at_y((*after_q).y())>(*after_q).x();
        if (q_larger_x == after_q_larger_x)
           return q_larger_x;
        else   // one smaller and one larger
        {
           // construct the other line
           Line_2 l_q = _construct_line_2(*q, *after_q); 
           if (l_q.is_horizontal())     // p is not horizontal
           {
              return (*q).x() > l_p.x_at_y((*q).y());
           }
           else 
           {
             bool p_larger_x = l_q.x_at_y((*p).y()) > (*p).x();
             bool after_p_larger_x = 
                  l_q.x_at_y((*after_p).y()) > (*after_p).x();

             CGAL_assertion (p_larger_x == after_p_larger_x);

             return !p_larger_x;
           }
        }
     }   
   }
}

}
