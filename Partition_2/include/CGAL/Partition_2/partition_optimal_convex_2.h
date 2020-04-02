// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>


// ===========================================================================
// There is a known bug in this algorithm that causes more than the optimal
// number of convex pieces to be reported in cases where there are many
// collinear vertices (as in a Hilbert polygon, for example).  More precisely,
// the problem is known to crop up in this situation:
//
//       5-----------4
//       |           |
//       |     2-----3
//       |     |
//       |     1-----0
//       |           |
//       |    13-----14
//       |     |
//       6    12--11
//       |        |
//       |     9--10
//       |     |
//       7-----8
//
// The problem arises because when decomposing the polygon from vertex 1 to 13
// point 2 is (quite correctly) indicated as not visible to point 13.  Thus
// it is believed an edge is necessary to divide the polygon (1 2 5 6 12 13).
//
// A hack that partially fixes this problem is implemented as follows:
// a vertex r is marked as visible from a point q for the purposes of
// the decompose function if p is the other endpoint of the edge containing r
// and p is visible from q.
//
// This causes the problem that decomposition from 8 to 12 indicates that
// valid vertices are 8 9 and 12. Diagonal (9 12) is valid, but (8 12) is
// also considered to be valid and necessary since 9 is a reflex vertex.
// The result is a polygon split with diagonals (2 5) (12 5) (9 12) and
// (13 1), which is obviously not optimal.
//
// To get around this problem, I currently postprocess during the cutting
// up of the polygon to remove unneeded diagonals and then
// achieve the optimal, but a better solution is surely desired....
//
// ===========================================================================

#ifndef CGAL_PARTITION_OPTIMAL_CONVEX_H
#define CGAL_PARTITION_OPTIMAL_CONVEX_H

#include <CGAL/license/Partition_2.h>


#include<CGAL/IO/Tee_for_output_iterator.h>
#include<CGAL/Partition_2/Partition_opt_cvx_edge.h>
#include<CGAL/Partition_2/Partition_opt_cvx_vertex.h>
#include<CGAL/Partition_2/Partition_opt_cvx_diagonal_list.h>
#include<CGAL/Partition_2/Matrix.h>
#include<CGAL/Partition_2/Turn_reverser.h>
#include<CGAL/Partition_2/Partitioned_polygon_2.h>
#include<CGAL/partition_is_valid_2.h>
#include<CGAL/Partition_traits_2.h>
#include<CGAL/Partition_2/partition_assertions.h>
#include<CGAL/Partition_2/Vertex_visibility_graph_2.h>
#include<utility>
#include<vector>
#include<iterator>
#include<iostream>


namespace CGAL {


#ifdef  CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
int partition_opt_cvx_debug_list_count = 0;
#endif


template <class Polygon, class Traits>
int partition_opt_cvx_best_so_far(Partition_opt_cvx_vertex& pivot_vertex,
                                  unsigned int extension,
                                  Polygon& polygon, const Traits& traits,
                                  Partition_opt_cvx_diagonal_list& diag_list)
{
#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
   std::cout << "best(" << pivot_vertex.vertex_num() << "."
             << partition_opt_cvx_debug_list_count
             << ", " << extension << ")" << std::endl;
#endif
   Partition_opt_cvx_stack_record best_so_far = pivot_vertex.best_so_far();
   while (!pivot_vertex.stack_empty()) {
      Partition_opt_cvx_stack_record old = pivot_vertex.stack_top();
#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
      std::cout << "best(" << pivot_vertex.vertex_num() << "."
                << partition_opt_cvx_debug_list_count << ", " << extension
                << ")" << " old = " << old.vertex_num()
                << ", " << old.value() << std::endl;
#endif
      typedef typename Traits::Left_turn_2    Left_turn_2;
      typedef typename Traits::Point_2       Point_2;
      Left_turn_2 left_turn = traits.left_turn_2_object();
      Turn_reverser<Point_2, Left_turn_2>  right_turn(left_turn);
      if (right_turn(polygon[old.vertex_num()],
                    polygon[pivot_vertex.vertex_num()], polygon[extension]))
      {
#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
         std::cout << "best(" << pivot_vertex.vertex_num() << "."
                   << partition_opt_cvx_debug_list_count
                   << ", " << extension << ")" << " returning(a) "
                   << best_so_far.value() << std::endl;
#endif
         diag_list = best_so_far.solution();
         return best_so_far.value();
      }
      else if (old.value() < best_so_far.value())
         best_so_far = old;
#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
      std::cout << "best(" << pivot_vertex.vertex_num() << "."
                << partition_opt_cvx_debug_list_count
                << ", " << extension << ") " << "popping off "
                << old.vertex_num() << ", " << old.value() << std::endl;
#endif
      pivot_vertex.stack_pop();
   }
#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
   std::cout << "best(" << pivot_vertex.vertex_num() << "."
             << partition_opt_cvx_debug_list_count
             << ", " << extension << ") returning(b) " << best_so_far.value()
             << std::endl;
#endif
   diag_list = best_so_far.solution();
#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
   std::cout << "     diagonal list " << diag_list << std::endl;
#endif
   return best_so_far.value();
}

template <class Polygon, class Traits>
void partition_opt_cvx_load(int current,
                            ::std::vector< Partition_opt_cvx_vertex >& v_list,
                            Polygon& polygon,
                            Matrix<Partition_opt_cvx_edge>& edges,
                            const Traits& traits)
{
    int previous;
    int num_polygons;
    Partition_opt_cvx_diagonal_list diag_list1, diag_list2;

#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
    std::cout << "load(" << v_list[current].vertex_num() << ")" << std::endl;
#endif
    for (previous = current-1; previous >= 0; previous--) {
#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
      std::cout << "load:  previous = " << v_list[previous].vertex_num()
                << std::endl;
#endif
      // must look at all valid edges and at all edges that are visible and
      // have something on the stack.  The latter check is necessary to make
      // sure solutions accumulate properly.
      if (edges[v_list[previous].vertex_num()]
               [v_list[current].vertex_num()].is_valid() ||
          (edges[v_list[previous].vertex_num()]
                [v_list[current].vertex_num()].is_visible() &&
           !v_list[previous].stack_empty()))
      {
         num_polygons = partition_opt_cvx_decompose(
                                  v_list[previous].vertex_num(),
                                  v_list[current].vertex_num(), polygon,
                                  edges, traits, diag_list1) +
                        partition_opt_cvx_best_so_far(v_list[previous],
                                                  v_list[current].vertex_num(),
                                                  polygon, traits, diag_list2);
         diag_list1.splice(diag_list1.end(), diag_list2);
#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
         std::cout << "load:  pushing previous = "
                   << v_list[previous].vertex_num() << " num_polygons = "
                   << num_polygons << " on stack "
                   << v_list[current].vertex_num() << "."
                   << partition_opt_cvx_debug_list_count << std::endl;
         std::cout << "     diagonal list = " << diag_list1 << std::endl;
#endif
         v_list[current].stack_push(v_list[previous].vertex_num(),
                                    num_polygons, diag_list1);
       }
    }
}

// pre: edge_num1 <= e_num <= edge_num2 but edge_num1 != edge_num2
template <class Polygon, class Traits>
bool collinearly_visible(unsigned int edge_num1, unsigned int e_num,
                         unsigned int edge_num2,
                         const Matrix<Partition_opt_cvx_edge>& edges,
                         const Polygon& polygon,
                         const Traits& traits)
{
   typedef typename Traits::Orientation_2                Orientation_2;
   typedef typename Traits::Point_2                      Point_2;
   Orientation_2 orientation = traits.orientation_2_object();

   if ((e_num == edge_num1+1 || e_num+1 == edge_num2) &&
       edges[edge_num1][edge_num2].is_visible() &&
       orientation(Point_2(polygon[edge_num1]), Point_2(polygon[e_num]),
                   Point_2(polygon[edge_num2])) ==
                  COLLINEAR)
     return true;
   else
     return false;
}


template <class Polygon, class Traits>
int partition_opt_cvx_decompose(unsigned int edge_num1, unsigned int edge_num2,
                                Polygon& polygon,
                                Matrix<Partition_opt_cvx_edge>& edges,
                                const Traits& traits,
                                Partition_opt_cvx_diagonal_list& diag_list)
{
#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
   std::cout << "decompose(" << edge_num1 << ", " << edge_num2 << ")";
#endif
   if (edges[edge_num1][edge_num2].is_done())  {
#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
      std::cout << " returning " << edges[edge_num1][edge_num2].value()
                << std::endl;
#endif
      diag_list = edges[edge_num1][edge_num2].solution();
      return edges[edge_num1][edge_num2].value();
   }
#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
   std::cout << std::endl;
#endif

   // temporarily invalidate this edge so we don't try to decompose on this
   // edge again
   Partition_opt_cvx_edge_validity old_validity;
   old_validity = edges[edge_num1][edge_num2].validity();
   edges[edge_num1][edge_num2].set_valid(PARTITION_OPT_CVX_NOT_VALID);

   std::vector< Partition_opt_cvx_vertex > v_list;
#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
   partition_opt_cvx_debug_list_count++;
#endif
   typedef typename Polygon::size_type  size_type;

   for (size_type e_num = edge_num1; e_num <= edge_num2; e_num++)
   {
       if ((edges[edge_num1][e_num].is_visible() &&
            edges[e_num][edge_num2].is_visible() ) ||
           collinearly_visible(edge_num1, static_cast<unsigned int>(e_num), edge_num2, edges, polygon,
                               traits) )
       {
         v_list.push_back(Partition_opt_cvx_vertex( static_cast<unsigned int>(e_num)));
       }
   }
   std::vector< int >::size_type v;

#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
   std::cout << "v_list(" << edge_num1 << ", " << edge_num2 << ")";
   for(v = 0; v < v_list.size(); v++) {
       std::cout << " " << v_list[v].vertex_num();
   }
   std::cout << std::endl;
#endif

   for(v = 0; v < v_list.size(); v++) {
       partition_opt_cvx_load(int(v), v_list, polygon, edges, traits);
   }

   int num_pieces = partition_opt_cvx_best_so_far(v_list[v_list.size()-1],
                                                  edge_num1, polygon, traits,
                                                  diag_list) + 1;
#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
   std::cout << "decompose: num_pieces = " << num_pieces << std::endl;
#endif
   edges[edge_num1][edge_num2].set_value(num_pieces);
   diag_list.push_back(Partition_opt_cvx_diagonal(edge_num1, edge_num2));
   edges[edge_num1][edge_num2].set_value(num_pieces);
   edges[edge_num1][edge_num2].set_solution(diag_list);
   edges[edge_num1][edge_num2].set_done(true);
   edges[edge_num1][edge_num2].set_valid(old_validity);
   // revalidate the edge; next time it will pick up the computed value
   // stored with this edge
#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
   std::cout << "decompose(" << edge_num1 << ", " << edge_num2 << "): "
             << " edge[" << edge_num1 << "][" << edge_num2 << "] set to "
             << edges[edge_num1][edge_num2] << std::endl;
   std::cout << " with diagonal list "
             << edges[edge_num1][edge_num2].solution()
             << std::endl;
   std::cout << "decompose(" << edge_num1 << ", " << edge_num2
             << "): returning " << num_pieces << std::endl;
   partition_opt_cvx_debug_list_count--;
#endif
   return num_pieces;
}

//
// implementation of the naive n^3 visibility algorithm
//
template <class Polygon, class Traits>
bool partition_opt_cvx_is_visible_n3(const Polygon& polygon, unsigned int i,
                                     unsigned int j, const Traits& traits)
{
   typedef typename Traits::Segment_2     Segment_2;
   typedef typename Polygon::size_type    size_type;
   typedef typename Traits::Left_turn_2    Left_turn_2;
   typedef typename Traits::Point_2       Point_2;
   typedef typename Traits::Construct_segment_2 Construct_segment_2;

   static Construct_segment_2 construct_segment_2 =
                                 traits.construct_segment_2_object();

   Segment_2 segment = construct_segment_2(polygon[i], polygon[j]);
   size_type prev_i = (i == 0)? polygon.size()-1: i - 1;
   size_type next_i = (i + 1)% polygon.size();
   size_type prev_j = (j == 0)? polygon.size()-1: j - 1;

   // determine if the edge goes through the interior of the polygon or not
   Left_turn_2 left_turn = traits.left_turn_2_object();
   Turn_reverser<Point_2, Left_turn_2> right_turn(left_turn);
   if (right_turn(polygon[prev_i], polygon[i], polygon[next_i]))
                                                      // concave angle
   {
     if (right_turn(polygon[prev_i], polygon[i], polygon[j]) &&
         right_turn(polygon[j], polygon[i], polygon[next_i]))
       return false;
   }
   else // left turn or straight
     if (right_turn(polygon[prev_i], polygon[i], polygon[j]) ||
         right_turn(polygon[j], polygon[i], polygon[next_i]))
       return false;

   size_type next_e;
   for (size_type e = 0; e < polygon.size(); e++) {
      if (e != i && e != prev_i && e != j && e != prev_j) {
         next_e = (e == polygon.size()-1)? 0 : e+1;
         Segment_2  edge = construct_segment_2(polygon[e], polygon[next_e]);
         if (do_intersect(segment, edge))
            return false;
      }
   }
   return true;
}

// when consecutive sequence of vertices are collinear, they must all be
// visible to each other as if there were no vertices in between.
template <class Polygon, class Traits>
void make_collinear_vertices_visible(Polygon& polygon,
                                     Matrix<Partition_opt_cvx_edge>& edges,
                                     const Traits& traits)
{
    typedef typename Polygon::size_type                   size_type;
    typedef typename Traits::Orientation_2                Orientation_2;
    typedef typename Traits::Point_2                      Point_2;
    Orientation_2 orientation = traits.orientation_2_object();

    size_type i;
    size_type prev_j, j;
    size_type k, next_k;

    // start at the beginning, move backwards as long as the points are
    // collinear; move forward as long as the points are collinear;
    // when you find the extremes make the larger one visible to the smaller
    // one loop until you reach the larger one each time starting again at the
    // larger
    i = polygon.size() - 1;
    prev_j = 0;
    j = 1;
    size_type start_i = 0;
    while (i > 0 &&
           orientation(Point_2(polygon[i]), Point_2(polygon[prev_j]), Point_2(polygon[j])) == COLLINEAR)
    {
       prev_j = i;
       start_i = i;
       i--;
    }
    i = 0;
    prev_j = 1;
    j = 2;
    while (j < polygon.size() &&
           orientation(Point_2(polygon[i]), Point_2(polygon[prev_j]), Point_2(polygon[j])) == COLLINEAR)
    {
       i++;
       prev_j++;
       j++;
    }
    // all points between start_i and prev_j are collinear so they must all
    // be visible to each other
    for (k = start_i; k != prev_j; )
    {
       next_k = k;
       do
       {
         next_k = (next_k == polygon.size() - 1) ? 0 : next_k+1;
         // the matrix should be upper triangular.
         if (k < next_k)
            edges[k][next_k].set_visible(true);
         else
            edges[next_k][k].set_visible(true);
       }
       while ( next_k != prev_j );
       k = (k == polygon.size() - 1) ? 0 : k+1;
    }
    i = prev_j;
    while (i < polygon.size())
    {
       prev_j = i+1;
       j = i+2;
       while (j < polygon.size() &&
              orientation(Point_2(polygon[i]), Point_2(polygon[prev_j]), Point_2(polygon[j])) ==
              COLLINEAR)
       {
           j++;
           prev_j++;
       }
       // the edge from the first collinear vertex to the last has
       // the same validity as the edge ending at the last collinear vertex
       if (prev_j < polygon.size())
          for (k = i; k != prev_j; k++)
          {
             next_k = k;
             do
             {
               next_k++;
               edges[k][next_k].set_visible(true);
             }
             while ( next_k != prev_j );
          }
       i = prev_j;
    }
}

template <class Polygon, class Traits>
void partition_opt_cvx_preprocessing(Polygon& polygon,
                                     Matrix<Partition_opt_cvx_edge>& edges,
                                     const Traits& traits)
{
    typedef typename Polygon::size_type                   size_type;

    typedef Vertex_visibility_graph_2<Traits>             Vis_graph;
    typedef typename Traits::Point_2                      Point_2;
    typedef std::pair<Point_2, Point_2>                   Point_pair;

    Vis_graph graph(polygon.begin(), polygon.end(), traits);

    size_type prev_i, i, next_i, next_next_i;
    size_type prev_j, j, next_j;

    for (i = 0; i < polygon.size(); i++)
    {
       prev_i = (i == 0)?polygon.size()-1:i-1;
       next_i = (i + 1)% polygon.size();
       next_next_i = (next_i + 1)% polygon.size();
       edges[i][i].set_visible(true);
       if (next_i != 0)
       {                                     // endpoints of edges are visible
          edges[i][next_i].set_visible(true);// and done (value == 0)
          edges[i][next_i].set_done(true);   // except for the last edge used
       }
       edges[i][next_i].set_valid(polygon[prev_i], polygon[i], polygon[next_i],
                             polygon[i], polygon[next_i], polygon[next_next_i],
                             traits);

       for (j = i + 2 ; j < polygon.size(); j++)
       {
          prev_j = j-1;
          if (graph.is_an_edge(Point_pair(polygon[i], polygon[j])))
          {
             next_j = (j + 1)% polygon.size();
             edges[i][j].set_visible(true);
             edges[i][j].set_valid(polygon[prev_i],polygon[i],polygon[next_i],
                                   polygon[prev_j],polygon[j],polygon[next_j],
                                   traits);
             if (j == i+2)
             {
                 edges[i][j].set_value(1);
                 Partition_opt_cvx_diagonal_list d;
                 d.push_back(Partition_opt_cvx_diagonal(static_cast<unsigned int>(i),
                                                        static_cast<unsigned int>(j)));
                 edges[i][j].set_solution(d);
                 edges[i][j].set_done(true);
             }
             // triangles are a base case.
          }
       }
    }
    make_collinear_vertices_visible(polygon, edges, traits);
}


template <class InputIterator, class OutputIterator, class Traits>
OutputIterator partition_optimal_convex_2(InputIterator first,
                                          InputIterator beyond,
                                          OutputIterator result,
                                          const Traits& traits)
{
   if (first == beyond) return result;

   typedef Partitioned_polygon_2< Traits >                 P_Polygon_2;
   typedef typename P_Polygon_2::iterator                  I;
   typedef typename P_Polygon_2::value_type                V;
   typedef typename P_Polygon_2::size_type                 S;
   typedef typename P_Polygon_2::difference_type           D;
   typedef Bidirectional_circulator_from_iterator<I,V,S,D> Circulator;


#if defined(CGAL_PARTITION_NO_POSTCONDITIONS) || \
    defined(CGAL_NO_POSTCONDITIONS)   || defined(NDEBUG)
   OutputIterator res(result);
#else
   typedef typename Traits::Polygon_2                      Polygon_2;
   Tee_for_output_iterator<OutputIterator, Polygon_2>      res(result);
#endif // no postconditions

   P_Polygon_2 polygon(first, beyond,traits);
   CGAL_partition_precondition(
    orientation_2(polygon.begin(), polygon.end(), traits) == COUNTERCLOCKWISE);

#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
   std::cout << "The polygon: " << std::endl;
   for (S i = 0; i < polygon.size(); i++)
      std::cout << polygon[i] << " ";
   std::cout << std::endl;
#endif

   Matrix<Partition_opt_cvx_edge> edges(polygon.size(), polygon.size());
   partition_opt_cvx_preprocessing(polygon, edges, traits);
#ifdef CGAL_PARTITION_OPTIMAL_CONVEX_DEBUG
   std::cout << "after preprocessing edges are (done, valid, visible, value): "
             << std::endl;
   std::cout << edges << std::endl;
#endif

   Partition_opt_cvx_diagonal_list diag_list;
   if (polygon.size() > 0)
   {
     partition_opt_cvx_decompose(0, static_cast<unsigned int>(polygon.size()-1), polygon, edges,
                                  traits, diag_list);

      diag_list.pop_back(); // the last diagonal added is the edge from last
                            // to first vertex (i.e., it is not a diagonal)
      Partition_opt_cvx_diagonal_list::const_iterator it;
      for (it = diag_list.begin(); it != diag_list.end(); it++)
      {
         Circulator source(polygon.begin(), polygon.end(),
                           polygon.begin() + (*it).first);
         Circulator target(polygon.begin(), polygon.end(),
                           polygon.begin() + (*it).second);
         polygon.insert_diagonal(source, target);
      }
      // the 1 is needed here to indicate that unnecessary edges should
      // be pruned away.  These crop up when there are collinear vertices.
      // See explanation at top of file.
      polygon.partition(res, 1);
      CGAL_partition_postcondition(
             convex_partition_is_valid_2(polygon.begin(), polygon.end(),
                                      res.output_so_far_begin(),
                                      res.output_so_far_end(), traits)
      );
   }

#if defined(CGAL_PARTITION_NO_POSTCONDITIONS) || \
    defined(CGAL_NO_POSTCONDITIONS)  || defined(NDEBUG)
   return res;
#else
   return res.to_output_iterator();
#endif // no postconditions

}

template <class InputIterator, class OutputIterator>
inline
OutputIterator partition_optimal_convex_2(InputIterator first,
                                          InputIterator beyond,
                                          OutputIterator result)
{
   typedef typename std::iterator_traits<InputIterator>::value_type Point_2;
   typedef typename Kernel_traits<Point_2>::Kernel  K;
   return partition_optimal_convex_2(first, beyond, result,
                                     Partition_traits_2<K>());
}

}

#endif // CGAL_PARTITION_OPTIMAL_CONVEX_H
