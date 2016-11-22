// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.5-I-117 $
// release_date  : $CGAL_Date: 2003/07/25 $
//
// file          : include/CGAL/partition_y_monotone_2.h
// package       : Partition_2 (1.53)
// maintainer    : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
// chapter       : Planar Polygon Partitioning
//
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: Y-monotone polygon partitioning
// ======================================================================

//
// Implementaion of the algorithm from pp 49--55 of "Computational Geometry 
// Algorithms and  Applications" by de Berg, van Kreveld, Overmars, and 
// Schwarzkopf for producing a partitioning of a polygon into y-monotone
// pieces.
//
// NOTE:  e_i = (v_i, v_{i+1})
//
// TREE:
//   "Therefore we store the edges of P intersecting the sweep line in the 
//    leaves of a dynamic binary search tree T.  The left-to-right order of
//    the leaves of T corresponds to the left-to-right order of the edges.
//    Because we are only interested in edges to the left of split and merge
//    vertices we only need to store edges in T that have the interior of P
//    to their right.  With each edge in T we store its helper."
//
//

#ifndef CGAL_PARTITION_Y_MONOTONE_H
#define CGAL_PARTITION_Y_MONOTONE_H

#include <CGAL/Indirect_not_less_yx_2.h>
#include <CGAL/Indirect_edge_compare.h>
#include <CGAL/IO/Tee_for_output_iterator.h>
#include <CGAL/partition_assertions.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Circulator_project.h>
#include <CGAL/function_objects.h>
#include <CGAL/functional.h>
#include <map>

#undef _DEBUG
#define _DEBUG 11
#include <CGAL/Nef_2/debug.h>

namespace CGAL {


enum Partition_y_mono_vertex_type {PARTITION_Y_MONO_START_VERTEX, 
                                   PARTITION_Y_MONO_SPLIT_VERTEX, 
                                   PARTITION_Y_MONO_REGULAR_VERTEX, 
                                   PARTITION_Y_MONO_COLLINEAR_VERTEX, 
                                   PARTITION_Y_MONO_MERGE_VERTEX, 
                                   PARTITION_Y_MONO_END_VERTEX};


//
// assumes CCW orientation of vertices
//
template <class BidirectionalCirculator, class Traits>
Partition_y_mono_vertex_type partition_y_mono_vertex_type(
                                BidirectionalCirculator c, 
                                const Traits& traits)
{
   BidirectionalCirculator previous = c;
   previous--;
   BidirectionalCirculator next = c;
   next++;
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   std::cout << "partition_y_mono__vertex_type: previous " << *previous 
             << " c " << *c << " next " << *next  << std::endl;
#endif
   typename Traits::Compare_y_2 compare_y_2 = traits.compare_y_2_object();
   typename Traits::Collinear_are_ordered_along_line_2 
     collinear_are_ordered_along_line_2
     = traits.collinear_are_ordered_along_line_2_object();

   if (compare_y_2(*previous, *c) == EQUAL &&
       compare_y_2(*next, *c) == EQUAL &&
       collinear_are_ordered_along_line_2( *previous, *c, *next))
      return PARTITION_Y_MONO_COLLINEAR_VERTEX;

   typename Traits::Less_yx_2   less_yx = traits.less_yx_2_object();
   typename Traits::Left_turn_2  left_turn = traits.leftturn_2_object();

   if (less_yx(*previous, *c)) 
   {
      if (less_yx(*next, *c))                // previous and next both less_yx
         if (left_turn(*previous, *c, *next)) // interior angle less than pi
             return PARTITION_Y_MONO_START_VERTEX;
         else                                // interior angle greater than pi
             return PARTITION_Y_MONO_SPLIT_VERTEX;
      else                                   // previous less_yx and next not
         return PARTITION_Y_MONO_REGULAR_VERTEX;
   }
   else 
   {
      if (less_yx(*c, *next))           // previous and next both not less_yx
        if (left_turn(*previous, *c, *next)) // interior angle less than pi
           return PARTITION_Y_MONO_END_VERTEX; 
        else                                // interior angle greater than pi
           return PARTITION_Y_MONO_MERGE_VERTEX;
      else                                 // next less_yx and previous not
        return PARTITION_Y_MONO_REGULAR_VERTEX;
   }
}

template <class Tree>
void partition_y_mono_print_tree(Tree tree)
{
   typedef typename Tree::iterator iterator;

   iterator it = tree.begin();
   for (; it != tree.end(); it++) {
    std::cout << "edge node " << *(*it).first << " helper " << *(*it).second 
              << std::endl;
   }
   std::cout << std::endl;
}

template <class BidirectionalCirculator, class Tree>
void partition_y_mono_handle_start_vertex(BidirectionalCirculator c, 
                                          Tree& tree)
{
   typedef typename Tree::value_type ValuePair;

#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   std::cout << *c << " is a start vertex " << std::endl;
#endif
   tree.insert(ValuePair(c, c));
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   std::cout << "partition_handle_start_vertex: after insert tree is " 
             << std::endl;
   partition_y_mono_print_tree(tree);
#endif
   // insert e_i (edge from *c to *++c) into "tree" with helper(e_i) = v_i
}

template <class BidirectionalCirculator, class Tree, 
          class OutputIterator, class Traits>
void partition_y_mono_handle_end_vertex(BidirectionalCirculator c, Tree& tree, 
                                        OutputIterator diagonals, 
                                        const Traits& traits )
{ 

#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   std::cout << *c << " is an end vertex " << std::endl;
#endif

   typedef typename Tree::iterator   tree_iterator;
   typedef typename Traits::Diagonal Diagonal;
   tree_iterator it;
   BidirectionalCirculator previous = c;
   previous--;

#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   std::cout << "partition_y_mono_handle_end_vertex: previous " << *previous 
             << std::endl;
#endif
   it = tree.find(previous);
   CGAL_assertion (it != tree.end());
   
   if (partition_y_mono_vertex_type((*it).second, traits) == 
          PARTITION_Y_MONO_MERGE_VERTEX) 
   {
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
       std::cout << "partition_y_mono_handle_end_vertex: diagonal " 
                 << *(*it).second << " to " << *c << std::endl;
#endif
       
       *diagonals++ = Diagonal( c, (*it).second);
   }
   tree.erase(it);
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   std::cout << "partition_y_mono_handle_end_vertex: after erase tree is " 
             << std::endl;
   partition_y_mono_print_tree(tree);
#endif
   // if helper(e_{i-1}) is a merge vertex
   //    insert diagonal connecting v_i to helper(e_{i-1})
   // delete e_{i-1} from tree
}

template<class BidirectionalCirculator, class Iterator, class Tree>
inline
void partition_y_mono_edge_directly_left(BidirectionalCirculator c, Tree& tree,
                                         Iterator& it)
{
   it = tree.lower_bound(c);  // edge directly to the left of v_i since the
                              // items in the tree are sorted from right to
                              // left
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   if (it != tree.end())
   std::cout << "partition_y_mono_edge_directly_left: lower_bound  edge node: "
             << *((*it).first) << " helper " << *((*it).second) << std::endl;
#endif
}

template <class BidirectionalCirculator, class Tree, class OutputIterator, class Traits>
void partition_y_mono_handle_split_vertex(BidirectionalCirculator c, 
					  Tree& tree, OutputIterator diagonals, 
					  const Traits& traits)
{
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   std::cout << *c << " is a split vertex " << std::endl;
#endif

   typedef typename Tree::iterator   tree_iterator;
   typedef typename Tree::value_type ValuePair;
   typedef typename Traits::Diagonal Diagonal;
   tree_iterator it;
   partition_y_mono_edge_directly_left(c, tree, it);
   if (it != tree.end()) 
   {
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
      std::cout << "partition_y_mono_handle_split_vertex: diagonal " 
                << *(*it).second << " to " << *c << std::endl;
#endif
      *diagonals++ = Diagonal( c, (*it).second);
      BidirectionalCirculator ej = (*it).first;
      tree.erase(it);
      tree.insert(ValuePair(ej, c));
   }
   tree.insert(ValuePair(c, c));
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   std::cout << "partition_y_mono_handle_split_vertex: "
             << "after erase and inserts tree is" << std::endl;
   partition_y_mono_print_tree(tree);
#endif
   // 1. find the edge e_j in tree directly to the left of v_i
   // 2. insert the diagonal connecting v_i to helper(e_j) 
   // 3. helper(e_j) = v_i
   // 4. Insert e_i in tree and set helper(e_i) to v_i
}

template <class BidirectionalCirculator, class Tree, 
          class OutputIterator, class Traits>
void partition_y_mono_handle_merge_vertex(BidirectionalCirculator c, 
                                          Tree& tree,
                                          OutputIterator diagonals, 
                                          const Traits& traits)
{
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   std::cout << *c << " is a merge vertex " << std::endl;
#endif

   typedef typename Tree::iterator     tree_iterator;
   typedef typename Tree::value_type   ValuePair;
   typedef typename Traits::Diagonal   Diagonal;
   BidirectionalCirculator prev = c;
   prev--;
   tree_iterator it = tree.find(prev);
   CGAL_assertion (it != tree.end());

   if (partition_y_mono_vertex_type((*it).second,traits) == 
         PARTITION_Y_MONO_MERGE_VERTEX)
   {
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
      std::cout << "partition_y_mono_handle_merge_vertex 1: diagonal " 
                << *(*it).second << " to " << *c << std::endl;
#endif
      *diagonals++ = Diagonal( c, (*it).second);
   }
   tree.erase(it);
   
   partition_y_mono_edge_directly_left(c, tree, it);
   if (it != tree.end()) 
   {
      if (partition_y_mono_vertex_type((*it).second,traits) == 
             PARTITION_Y_MONO_MERGE_VERTEX) 
      {
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
      std::cout << "partition_y_mono_handle_merge_vertex 2: diagonal " 
                << *(*it).second << " to " << *c << std::endl;
#endif
      *diagonals++ = Diagonal( c, (*it).second);
      }
      BidirectionalCirculator ej = (*it).first;
      tree.erase(it);
      tree.insert(ValuePair(ej,c));
   }
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   std::cout << "partition_y_mono_handle_merge_vertex: after erase(s) tree is "
             << std::endl;
   partition_y_mono_print_tree(tree);
#endif
   // 1. if helper(e_{i-1}) is a merge vertex
   //       insert the diagonal connecting v_i to helper(e_{i-1})
   // 2.  delete e_{i-1} from tree
   // 3.  find the edge e_j in tree directly to the left of v_i
   // 4.  if helper(e_j) is a merge vertex
   //        insert diagonal connecting v_i to helper(e_j) in polygon
   // 5.  helper(e_j) = v_i
}

template <class BidirectionalCirculator, class Traits>
bool partition_y_mono_interior_to_right(BidirectionalCirculator c,
                                        const Traits& traits)
{
   typename Traits::Compare_y_2 compare_y_2 = traits.compare_y_2_object();

   BidirectionalCirculator previous = c; previous--;

   Comparison_result cmp_y = compare_y_2(*previous, *c);
   if (cmp_y == LARGER) return true;

   BidirectionalCirculator next = c; next++;

   if (cmp_y == EQUAL && compare_y_2(*next, *c) == SMALLER) return true;

   return false;
}

template <class BidirectionalCirculator, class Tree, class OutputIterator,
          class Traits>
void partition_y_mono_handle_regular_vertex(BidirectionalCirculator c, 
                                            Tree& tree, 
                                            OutputIterator diagonals, 
                                            const Traits& traits )
{
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   std::cout << *c << " is a regular vertex " << std::endl;
#endif
   typedef typename Tree::iterator     tree_iterator;
   typedef typename Tree::value_type   ValuePair;
   typedef typename Traits::Diagonal   Diagonal;
   tree_iterator it;

   BidirectionalCirculator previous = c;
   previous--;

   if (partition_y_mono_interior_to_right(c, traits))  
   {
      it = tree.find(previous);
      CGAL_assertion( it != tree.end() );

      if (partition_y_mono_vertex_type((*it).second, traits) == 
             PARTITION_Y_MONO_MERGE_VERTEX)
      {
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
         std::cout << "partition_y_mono_handle_regular_vertex 1: diagonal " 
                   << *(*it).second << " to " << *c << std::endl;
#endif
         *diagonals++ = Diagonal( c, (*it).second);
      }
      tree.erase(it);
      tree.insert(ValuePair(c,c));
   }
   else 
   {
      partition_y_mono_edge_directly_left(c, tree, it);
      CGAL_assertion (it != tree.end());

      if (partition_y_mono_vertex_type((*it).second, traits) == 
             PARTITION_Y_MONO_MERGE_VERTEX)
      {
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
         std::cout << "partition_y_mono_handle_regular_vertex 2: diagonal " 
                   << *c << " to " << *(*it).second << std::endl;
#endif
         *diagonals++ = Diagonal( c, (*it).second);
      }
      BidirectionalCirculator ej = (*it).first;
      tree.erase(it);
      tree.insert(ValuePair(ej,c));
   }
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   std::cout << "partition_y_mono_handle_regular_vertex: "
             << "after erase and insert tree is" << std::endl;
   partition_y_mono_print_tree(tree);
#endif
   //  if interior of polygon lies to the right of v_i
   //     if helper(e_{i-1}) is a merge vertex
   //        insert diagonal connecting v_i to helper(e_{i-1}) in polygon
   //     delete e_{i-1} from tree
   //     insert e_i in tree and set helper(e_i) to v_i
   //  else
   //     find the edge e_j in tree directly left of v_i
   //     if helper(e_j) is a merge vertex
   //        insert diagonal connecting v_i to helper(e_j) in D
   //     helper(e_j) = v_i
}

template <class BidirectionalCirculator, class Tree>
void partition_y_mono_handle_collinear_vertex(BidirectionalCirculator c, 
                                              Tree& tree)
{
   typedef typename Tree::iterator     tree_iterator;
   typedef typename Tree::value_type   ValuePair;
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   std::cout << *c << " is a collinear vertex " << std::endl;
#endif

   tree_iterator it;

   BidirectionalCirculator previous = c;
   previous--;
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   std::cout << *previous << " is the previous vertex " << std::endl;
#endif

   it = tree.find(previous);
   if ( it != tree.end() )
   {
#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
      std::cout << "partition_y_mono_handle_collinear_vertex : removing " 
                << *(*it).first << std::endl;
#endif
      tree.erase(it);
   }
   tree.insert(ValuePair(c,c));
}


template <class InputIterator, class OutputIterator, class Traits>
void  partition_y_monotone_2( InputIterator first, InputIterator beyond, 
			      OutputIterator diagonals, const Traits& traits)
{
   typedef typename InputIterator::value_type Circulator;
#ifdef _DEBUG
   typedef Container_from_circulator<Circulator> Container;
   typedef typename Traits::Point_2 Point_2;
#endif

   std::vector<Circulator> circulators;
   InputIterator bi = first;
   Circulator circ(*bi), done(circ);
#ifdef _DEBUG
   //CGAL_assertion_code(Container cycle(circ));
   //CGAL_NEF_TRACEN("partition_y_monotone_2(): outer cycle:");
   //std::copy( cycle.begin(), cycle.end(), std::ostream_iterator<Point_2>( std::cerr, "\n"));
#endif
   //CGAL_partition_precondition
   //CGAL_partition_warning
   //  (orientation_2( cycle.begin(), cycle.end(), traits) == COUNTERCLOCKWISE);
   CGAL_For_all( circ, done){
     circulators.push_back(circ);
   }
   for( ++bi; bi != beyond; ++bi) {
     Circulator circ(*bi), done(circ);
#ifdef _DEBUG
     //CGAL_assertion_code(Container cycle(circ));
     //CGAL_NEF_TRACEN("partition_y_monotone_2(): inner cycle:");
     //std::copy( cycle.begin(), cycle.end(), std::ostream_iterator<Point_2>( std::cerr, "\n"));
#endif
     //CGAL_partition_precondition
     //CGAL_partition_warning
     //  (orientation_2( cycle.begin(), cycle.end(), traits) == CLOCKWISE);
     CGAL_For_all( circ, done){
       circulators.push_back(circ);
     }
   }

   std::sort(circulators.begin(), circulators.end(), 
	     Indirect_not_less_yx_2<Traits>(traits));

#ifdef CGAL_PARTITION_Y_MONOTONE_DEBUG
   std::cout << "Initial vertex list: ";
   for( typename std::vector<Circulator>::const_iterator it = circulators.begin();
       it != circulators.end();
       it++){
     std::cout << **it << " " ; 
   }
   std::cout << std::endl;
#endif

   typedef std::map<Circulator, Circulator, 
                    Indirect_edge_compare<Circulator, Traits> > Tree;
   Tree tree;

   typename std::vector<Circulator>::iterator it = circulators.begin();
   for (; it != circulators.end(); it++) {
      switch (partition_y_mono_vertex_type(*it, traits)) 
      {
         case PARTITION_Y_MONO_START_VERTEX:
            partition_y_mono_handle_start_vertex(*it, tree);
            break;
         case PARTITION_Y_MONO_SPLIT_VERTEX:
            partition_y_mono_handle_split_vertex(*it, tree, diagonals, traits);
            break;
         case PARTITION_Y_MONO_END_VERTEX:
            partition_y_mono_handle_end_vertex(*it, tree, diagonals, traits);
            break;
         case PARTITION_Y_MONO_MERGE_VERTEX:
            partition_y_mono_handle_merge_vertex(*it, tree, diagonals, traits);
            break;
         case PARTITION_Y_MONO_REGULAR_VERTEX:
            partition_y_mono_handle_regular_vertex(*it, tree, diagonals, traits);
            break;
         case PARTITION_Y_MONO_COLLINEAR_VERTEX:
            partition_y_mono_handle_collinear_vertex(*it, tree);
            break;
      }
   }

   /*
   CGAL_partition_postcondition(
       y_monotone_partition_is_valid_2(polygon.begin(), polygon.end(),
                                       res.output_so_far_begin(),
                                       res.output_so_far_end(), traits));
   */
}

template <class InputIterator, class OutputIterator>
inline void partition_y_monotone_2( InputIterator first, 
				    InputIterator beyond, 
				    OutputIterator diagonals)
{
  typedef typename std::iterator_traits<InputIterator>::value_type 
    BoundaryCirculator;
  typedef typename std::iterator_traits<BoundaryCirculator>::value_type
    Point_2;
   typedef typename Kernel_traits<Point_2>::Kernel K;

   partition_y_monotone_2( first, beyond, diagonals, 
			   Partition_traits_2< BoundaryCirculator, K>());
}

}

#endif // CGAL_PARTITION_Y_MONOTONE_H
