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
// file          : include/CGAL/Partitioned_polygon_2.C
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
// implementation: Polygon with diagonals used to represent a partitioning
// ============================================================================

#include <CGAL/Partitioned_polygon_2.h>

namespace CGAL {

template <class Traits>
void 
Partitioned_polygon_2<Traits>::prune_diagonals()
{
    Circulator first(begin(), end(), begin());
    Circulator c = first;
    typedef Partition_vertex<Traits>::Diagonal_iterator Diagonal_iterator;

    Diagonal_iterator d;
#ifdef CGAL_PARTITIONED_POLY_DEBUG
    std::cout << "pruning diagonals ..." << std::endl;
#endif
    do {
       d = (*c).diagonals_begin();
       while (d != (*c).diagonals_end()) {
          if (!diagonal_is_necessary(c, *d)) 
          {
#ifdef CGAL_PARTITIONED_POLY_DEBUG
             std::cout << "   removing from " << *c << " to " << **d 
                       << std::endl;
#endif
             (**d).diagonal_erase(c);
             d = (*c).diagonal_erase(d);
          }
          else
          {
            d++;
          }
       }
#ifdef CGAL_PARTITIONED_POLY_DEBUG
       (*c).print_diagonals();
#endif
       (*c).reset_current_diagonal();  
    }
    while (++c != first);
}

// the pruning is probably no longer necessary
template <class Traits>
template <class OutputIterator>
OutputIterator 
Partitioned_polygon_2<Traits>::partition(OutputIterator result, bool prune) 
{
    // walk through each vertex and sort the diagonals 
    Circulator first(begin(), end());
    Circulator c = first;
    Circulator next;
    Circulator prev = c;
    prev--;
    do 
    {
       next = c;
       next++;
       (*c).sort_diagonals(prev, next);
#ifdef CGAL_PARTITIONED_POLY_DEBUG
       (*c).print_diagonals();
#endif
       prev = c;
    }
    while (++c != first);

    // now remove any diagonals that do not cut a reflex angle at one end
    if (prune) prune_diagonals();

#ifdef CGAL_PARTITIONED_POLY_DEBUG
    c = first;
    do 
    {
       (*c).print_diagonals();
    }
    while (++c != first);
#endif

    make_polygon(first, result);
    return result;
}

template <class Traits>
template <class OutputIterator>
Partitioned_polygon_2<Traits>::Circulator
Partitioned_polygon_2<Traits>::make_polygon(
                             Partitioned_polygon_2<Traits>::Circulator start, 
                             OutputIterator& result)
{
    Subpolygon_2 new_polygon;
    Circulator next = start;
    do 
    {
       new_polygon.push_back(*next);
#ifdef CGAL_PARTITIONED_POLY_DEBUG
       std::cout << "adding vertex " << *next << std::endl;
#endif
       Circulator diag;
       if ((*next).has_unused_diagonals()) 
       {
          diag = (*next).current_diagonal();
#ifdef CGAL_PARTITIONED_POLY_DEBUG
          std::cout << "diagonal endpoint: " << *diag << std::endl;
#endif
          (*next).advance_diagonal();
          if (diag == start) 
          {
             *result = new_polygon;
             result++;
             return next;
          }
          else 
          {
             next = make_polygon(next, result);
          }
       }
       else next++;
    } while (next != start);
    *result = new_polygon;
    result++;
    return next;
    // if there are no diagonals at this vertex
    //    push on the vertex
    // else if the first diagonal closes the polygon
    //    close the polygon
    //    return the current vertex (NOT the other end of the diagonal)
    // else 
    //    remove the first diagonal
    //    recur, starting a new polygon at this vertex and return the 
    //      vertex where the new polygon ended
    //    continue from the last vertex of the new polygon
}

template <class Traits>
bool 
Partitioned_polygon_2<Traits>::cuts_reflex_angle(
            Partitioned_polygon_2<Traits>::Circulator vertex_ref, 
            Partitioned_polygon_2<Traits>::Circulator diag_endpoint)
{
   Circulator prev = vertex_ref; prev--;
   Circulator next = vertex_ref; next++;

   typedef Partition_vertex<Traits>::Diagonal_iterator Diagonal_iterator;
   
   // find diag_endpoint in vertex_ref's list of diagonals
   Diagonal_iterator d_it;
   for (d_it = (*vertex_ref).diagonals_begin(); 
        d_it != (*vertex_ref).diagonals_end() && diag_endpoint != *d_it; 
        d_it++)
   {
      prev = *d_it;
   }
   Diagonal_iterator next_d_it = d_it;
   next_d_it++;
   if (next_d_it == (*vertex_ref).diagonals_end()) 
   {
      next = vertex_ref;
      next++;
   }
   else
      next = *next_d_it;

   return _rightturn(*prev, *vertex_ref, *next);
}

}
