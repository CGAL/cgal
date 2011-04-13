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
// file          : include/CGAL/partition_optimal_convex_2.h
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
// implementation: Optimal convex polygon partitition
// ============================================================================

#ifndef CGAL_PARTITION_OPTIMAL_CONVEX_H
#define CGAL_PARTITION_OPTIMAL_CONVEX_H

#include<CGAL/IO/Tee_for_output_iterator.h>
#include<CGAL/Partition_opt_cvx_edge.h>
#include<CGAL/Partition_opt_cvx_vertex.h>
#include<CGAL/Partition_opt_cvx_diagonal_list.h>
#include<CGAL/Matrix.h>
#include<CGAL/Turn_reverser.h>
#include<CGAL/Partitioned_polygon_2.h>
#include<CGAL/partition_is_valid_2.h>
#include<CGAL/Partition_traits_2.h>
#include<CGAL/partition_assertions.h>
#include<CGAL/Vertex_visibility_traits_2.h>
#include<CGAL/Vertex_visibility_graph_2.h>
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
      typedef typename Traits::Leftturn_2    Leftturn_2;
      typedef typename Traits::Point_2       Point_2;
      Leftturn_2 leftturn = traits.leftturn_2_object();
      Turn_reverser<Point_2, Leftturn_2>  rightturn(leftturn);
      if (rightturn(polygon[old.vertex_num()], 
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
       if (edges[edge_num1][e_num].is_visible() && 
           edges[e_num][edge_num2].is_visible()) 
       {
          v_list.push_back(Partition_opt_cvx_vertex(e_num));
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
   typedef typename Traits::R             R;
   typedef typename R::Segment_2          Segment_2;
   typedef typename Polygon::size_type    size_type;
   typedef typename Traits::Leftturn_2    Leftturn_2;
   typedef typename Traits::Point_2       Point_2;
   typedef typename Traits::Construct_segment_2 Construct_segment_2;

   static Construct_segment_2 construct_segment_2 = 
                                 traits.construct_segment_2_object();

   Segment_2 segment = construct_segment_2(polygon[i], polygon[j]);
   size_type prev_i = (i == 0)? polygon.size()-1: i - 1;
   size_type next_i = (i + 1)% polygon.size();
   size_type prev_j = (j == 0)? polygon.size()-1: j - 1;

   // determine if the edge goes through the interior of the polygon or not
   Leftturn_2 leftturn = traits.leftturn_2_object();
   Turn_reverser<Point_2, Leftturn_2> rightturn(leftturn);
   if (rightturn(polygon[prev_i], polygon[i], polygon[next_i]))// concave angle
   {
     if (rightturn(polygon[prev_i], polygon[i], polygon[j]) &&
         rightturn(polygon[j], polygon[i], polygon[next_i]))
       return false;
   }
   else // left turn or straight
     if (rightturn(polygon[prev_i], polygon[i], polygon[j]) ||
         rightturn(polygon[j], polygon[i], polygon[next_i]))
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


template <class Polygon, class Traits>
void partition_opt_cvx_preprocessing(Polygon& polygon, 
                                     Matrix<Partition_opt_cvx_edge>& edges, 
                                     const Traits& traits)
{
    typedef typename Polygon::size_type                   size_type;

    typedef typename Traits::R                            R;
    typedef typename Polygon::iterator                    Vertex_iterator;
    typedef Vertex_visibility_traits_2<R>                 Vis_traits;
    typedef Vertex_visibility_graph_2<Traits>             Vis_graph;
    typedef typename Traits::Point_2                      Point_2;
    typedef std::pair<Point_2, Point_2>                   Point_pair;

    Vis_graph graph(polygon.begin(), polygon.end());

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
          if (graph.is_an_edge(Point_pair(polygon[i], polygon[j])))
          {
             prev_j = (j == 0)?polygon.size()-1:j-1;
             next_j = (j + 1)% polygon.size();
             edges[i][j].set_visible(true);
             edges[i][j].set_valid(polygon[prev_i],polygon[i],polygon[next_i],
                                   polygon[prev_j],polygon[j],polygon[next_j],
                                   traits);
             if (j == i+2) 
             {
                 edges[i][j].set_value(1); 
                 Partition_opt_cvx_diagonal_list d;
                 d.push_back(Partition_opt_cvx_diagonal(i,j));
                 edges[i][j].set_solution(d); 
                 edges[i][j].set_done(true); 
             }
             // triangles are a base case.
          }
       }
   }
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

   P_Polygon_2 polygon(first, beyond);
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
      partition_opt_cvx_decompose(0, polygon.size()-1, polygon, edges, 
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
      polygon.partition(res, 0);
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
   return CGAL_partition_optimal_convex_2(first, beyond, result, 
                                          reinterpret_cast<Point_2*>(0));
}

template <class InputIterator, class OutputIterator, class R>
inline
OutputIterator CGAL_partition_optimal_convex_2(InputIterator first, 
                                               InputIterator beyond,
                                               OutputIterator result, 
                                               Point_2<R>*)
{
   return partition_optimal_convex_2(first, beyond, result,
                                     Partition_traits_2<R>());
}


}

#endif // CGAL_PARTITION_OPTIMAL_CONVEX_H
