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
// file          : include/CGAL/partition_is_valid_2.h
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
// implementation: Polygon partitioning validity checking functions.
// ============================================================================

#ifndef CGAL_PARTITION_IS_VALID_2_H
#define CGAL_PARTITION_IS_VALID_2_H

#include <list>
#include <utility>
#include <iterator>
#include <CGAL/partition_assertions.h>
#include <CGAL/Partitioned_polygon_2.h>
#include <CGAL/Partition_vertex_map.h>
#include <CGAL/ch_selected_extreme_points_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Partition_is_valid_traits_2.h>
#include <CGAL/Polygon_2.h>

// NOTE:  this could possibly be checked using a planar map overlay, but
//        then the traits class would have to require lots of other things
//        and you have to do many overlaps, not just one so it is
//        less efficient.

namespace CGAL {

template <class InputIterator, class Vertex_list>
bool 
polygons_equal(InputIterator first, InputIterator last, Vertex_list vertices)
{
   typedef typename Vertex_list::iterator        I;
   typedef Circulator_from_iterator<I>  Circulator;

   Circulator v_first(vertices.begin(), vertices.end(), vertices.begin());
   Circulator v_circ = v_first;

   // find the first (input) vertex in the list of vertices
   for (v_circ = v_first; *v_circ != *first && ++v_circ != v_first;) {}

   // now look through both lists until you find a vertex that is not
   // the same or you reach the end of the vertices
   while (first != last && *v_circ == *first)
   {
#ifdef CGAL_PARTITION_CHECK_DEBUG
      std::cout << *first << " is in the right place " << std::endl;
#endif
      v_circ++; first++;
   }

   if (first == last) return true;
#ifdef CGAL_PARTITION_CHECK_DEBUG
   std::cout << *first << " is out of order. ";
#endif
   return false;
}


template<class InputIterator1, class InputIterator2, class Traits> 
bool
partition_is_valid_2 (InputIterator1 point_first, InputIterator1 point_last,
                      InputIterator2 poly_first, InputIterator2 poly_last,
                      const Traits& traits) 
{
   if (poly_first == poly_last)  return (point_first == point_last);

   typedef typename Traits::Polygon_2::Vertex_const_iterator   
                                                       Poly_vtx_const_iterator;
   typedef typename Traits::Point_2                    Point_2;
   typedef Edge_list<Traits>                           Edge_list;
   typedef std::pair<Point_2, Edge_list>               P_Vertex;
   typedef Partition_vertex_map<Traits>                P_Vertex_map;
   typedef typename P_Vertex_map::iterator             Set_iterator;
   typedef std::pair<Set_iterator, bool>               Set_iterator_pair;
   typedef typename Traits::Is_valid                   Is_valid;

   Poly_vtx_const_iterator begin, end, v_it;

   Is_valid is_valid = traits.is_valid_object(traits);

   std::list<Point_2>    orig_poly;
   for (;point_first != point_last; point_first++)
      orig_poly.push_back(*point_first);

   CGAL_partition_precondition(orientation_2(orig_poly.begin(),orig_poly.end(),
                                             traits) == COUNTERCLOCKWISE);
   P_Vertex_map  output_vertex_set;
   Set_iterator_pair v_loc_pair, begin_v_loc_pair, prev_v_loc_pair;

   int poly_num = 0;
   for (; poly_first != poly_last; poly_first++, poly_num++)
   {
      begin = (*poly_first).vertices_begin();
      end = (*poly_first).vertices_end();
#ifdef CGAL_PARTITION_CHECK_DEBUG
         std::cout << "Polygon " << poly_num << " is " << std::endl;
         std::cout << *poly_first << std::endl;
#endif
      CGAL_partition_assertion (
           orientation_2(begin, end, traits) == COUNTERCLOCKWISE);
      if (!is_valid(begin, end)) 
      {
#ifdef CGAL_PARTITION_CHECK_DEBUG
         std::cout << "It does NOT have the tested property." << std::endl;
#endif
         return false;
      }
      begin_v_loc_pair=output_vertex_set.insert(P_Vertex(*begin, Edge_list()));
      prev_v_loc_pair = begin_v_loc_pair;
      v_it = begin;
      for (v_it++; v_it != end; v_it++)
      {
         v_loc_pair = output_vertex_set.insert(P_Vertex(*v_it, Edge_list()));
         output_vertex_set.insert_next_edge(prev_v_loc_pair.first, 
                                            v_loc_pair.first,
                                            poly_num);
         output_vertex_set.insert_prev_edge(v_loc_pair.first, 
                                            prev_v_loc_pair.first,
                                            poly_num);
         prev_v_loc_pair = v_loc_pair;
      }
      output_vertex_set.insert_next_edge(prev_v_loc_pair.first, 
                                         begin_v_loc_pair.first,
                                         poly_num);
      output_vertex_set.insert_prev_edge(begin_v_loc_pair.first, 
                                         prev_v_loc_pair.first,
                                         poly_num);
   }
   if (output_vertex_set.polygons_overlap()) return false;

   std::list<Point_2>  union_polygon;
   output_vertex_set.union_vertices(std::back_inserter(union_polygon));

#ifdef CGAL_PARTITION_CHECK_DEBUG
   std::list<Point_2>::iterator  poly_iterator;
   std::cout << "union polygon is " << std::endl;
   for (poly_iterator = union_polygon.begin(); 
        poly_iterator != union_polygon.end(); poly_iterator++)
   {
      std::cout << *poly_iterator << " ";
   }
   std::cout << std::endl;
#endif

   if (!polygons_equal(orig_poly.begin(), orig_poly.end(), union_polygon)) 
      return false;

   return true;
}

template<class InputIterator1, class InputIterator2>
bool
partition_is_valid_2 (InputIterator1 point_first, InputIterator1 point_last,
                      InputIterator2 poly_first, InputIterator2 poly_last)
{
   typedef typename std::iterator_traits<InputIterator1>::value_type   Point_2;
   return CGAL_partition_is_valid_2(point_first, point_last,
                                    poly_first, poly_last,
                                    reinterpret_cast<Point_2*>(0));
}

template<class InputIterator1, class InputIterator2, class R>
bool
CGAL_partition_is_valid_2 (InputIterator1 point_first, 
                           InputIterator1 point_last,
                           InputIterator2 poly_first, 
                           InputIterator2 poly_last,
                           Point_2<R>*)
{
   typedef Partition_traits_2<R>                   Traits;
   typedef Is_vacuously_valid<Traits>              Is_valid;

   Partition_is_valid_traits_2<Traits, Is_valid>   validity_traits;

   return partition_is_valid_2(point_first, point_last, poly_first, poly_last,
                               validity_traits);
}


template<class InputIterator1, class InputIterator2, class Traits>
bool 
convex_partition_is_valid_2(InputIterator1 point_first,
                            InputIterator1 point_last,
                            InputIterator2 poly_first,
                            InputIterator2 poly_last,
                            const Traits& traits)
{
   typedef typename Traits::Is_convex_2                 Is_convex_2;
   Partition_is_valid_traits_2<Traits, Is_convex_2>     validity_traits;

   return partition_is_valid_2(point_first, point_last, poly_first, poly_last,
                               validity_traits);
}

template<class InputIterator1, class InputIterator2>
bool 
convex_partition_is_valid_2(InputIterator1 point_first,
                            InputIterator1 point_last,
                            InputIterator2 poly_first,
                            InputIterator2 poly_last)
{
   typedef typename std::iterator_traits<InputIterator1>::value_type   Point_2;
   return CGAL_convex_partition_is_valid_2(point_first, point_last, 
                                           poly_first, poly_last, 
                                           reinterpret_cast<Point_2*>(0));
}

template<class InputIterator1, class InputIterator2, class R>
bool 
CGAL_convex_partition_is_valid_2(InputIterator1 point_first,
                                 InputIterator1 point_last,
                                 InputIterator2 poly_first,
                                 InputIterator2 poly_last,
                                 Point_2<R>*)
{
   return convex_partition_is_valid_2(point_first, point_last, 
                                      poly_first, poly_last, 
                                      Partition_traits_2<R>());
}

template<class InputIterator1, class InputIterator2, class Traits>
bool 
y_monotone_partition_is_valid_2(InputIterator1 point_first,
                                InputIterator1 point_last,
                                InputIterator2 poly_first,
                                InputIterator2 poly_last,
                                const Traits& traits)
{
   typedef typename Traits::Is_y_monotone_2                Is_y_monotone_2;

   Partition_is_valid_traits_2<Traits, Is_y_monotone_2>    validity_traits;

   return partition_is_valid_2(point_first, point_last, poly_first, poly_last,
                               validity_traits);
}

template<class InputIterator1, class InputIterator2>
bool 
y_monotone_partition_is_valid_2(InputIterator1 point_first,
                                InputIterator1 point_last,
                                InputIterator2 poly_first,
                                InputIterator2 poly_last)
{
   typedef typename std::iterator_traits<InputIterator1>::value_type   Point_2;
   return CGAL_y_monotone_partition_is_valid_2(point_first, point_last, 
                                               poly_first, poly_last,
                                               reinterpret_cast<Point_2*>(0));
}

template<class InputIterator1, class InputIterator2, class R>
bool
CGAL_y_monotone_partition_is_valid_2(InputIterator1 point_first,
                                     InputIterator1 point_last,
                                     InputIterator2 poly_first,
                                     InputIterator2 poly_last,
                                     Point_2<R>*)
{
   return y_monotone_partition_is_valid_2(point_first, point_last, 
                                          poly_first, poly_last, 
                                          Partition_traits_2<R>());
}

}

#endif // CGAL_PARTITION_IS_VALID_2_H
