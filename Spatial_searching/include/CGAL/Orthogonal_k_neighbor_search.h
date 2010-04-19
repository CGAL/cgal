// Copyright (c) 2002 Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef CGAL_ORTHOGONAL_K_NEIGHBOR_SEARCH_H
#define CGAL_ORTHOGONAL_K_NEIGHBOR_SEARCH_H

#include <cstring>
#include <list>
#include <queue>
#include <set>
#include <memory>
#include <CGAL/Kd_tree_node.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Splitters.h>

namespace CGAL {

template <class SearchTraits, 
          class Distance_= Euclidean_distance<SearchTraits>,
          class Splitter_= Sliding_midpoint<SearchTraits> ,
          class Tree_= Kd_tree<SearchTraits, Splitter_, Tag_true> >
class Orthogonal_k_neighbor_search {

public:

  typedef Splitter_ Splitter;
  typedef Tree_  Tree;
  typedef Distance_ Distance;
  typedef typename SearchTraits::Point_d Point_d;
  typedef typename Distance::Query_item Query_item;

  typedef typename SearchTraits::FT FT;
  typedef std::pair<Point_d,FT> Point_with_transformed_distance;

  typedef typename Tree::Node_handle Node_handle;

  typedef typename Tree::Point_d_iterator Point_d_iterator;

private:

  // Comparison functor to sort a set of points
  // in increasing or decreasing order (key is distance).
  class Distance_larger 
  {
    bool search_nearest;

  public:

    Distance_larger(bool search_the_nearest_neighbour) 
      : search_nearest(search_the_nearest_neighbour)
    {}

    bool operator()(const Point_with_transformed_distance& p1, 
                    const Point_with_transformed_distance& p2) const  
    {
      if (search_nearest)
        return (p1.second < p2.second);
      else
        return (p2.second < p1.second);
    }
  };
  
  // Set of points, sorted by distance, in increasing or decreasing order.
  typedef std::multiset<Point_with_transformed_distance, Distance_larger> NN_list;

public:

  typedef typename NN_list::const_iterator iterator;

private:

  int number_of_internal_nodes_visited;
  int number_of_leaf_nodes_visited;
  int number_of_items_visited;

  bool search_nearest;

  FT multiplication_factor;
  Query_item query_object;
  int total_item_number;
  FT distance_to_root;   

  NN_list l; // Set of points, sorted by distance
  unsigned int max_k;
  
  Distance distance_instance;

private:

  // Test if we should continue searching 
  inline bool branch(FT distance) 
  {
    if (l.size()<max_k) 
      return true;
    else 
      if (search_nearest) 
        return (distance*multiplication_factor < (--l.end())->second);
      else 
        return (distance > l.begin()->second*multiplication_factor);
  }

  // Try to insert point *I
  void insert(Point_d* I, FT dist) 
  {
    // Shall we insert I?
    bool insert;
    if (l.size()<max_k) 
      insert=true;
    else if (search_nearest) 
      insert = ( dist < (--l.end())->second ); 
    else 
      insert=(dist > (--l.end())->second);
      
    if (insert) 
    {
      Point_with_transformed_distance NN_Candidate(*I,dist);
      l.insert(NN_Candidate);
      if (l.size() > max_k)
        l.erase(--l.end()); 
    }
  }


public:

  iterator begin() const
  {
    return l.begin();
  }

  iterator end() const
  {
    return l.end();
  }


  // constructor
  Orthogonal_k_neighbor_search(Tree& tree, const Query_item& q,  
    unsigned int k=1, FT Eps=FT(0.0), bool Search_nearest=true, const Distance& d=Distance())
    : number_of_internal_nodes_visited(0), number_of_leaf_nodes_visited(0), number_of_items_visited(0), 
    search_nearest(Search_nearest), multiplication_factor(d.transformed_distance(1.0+Eps)), query_object(q), 
    total_item_number(tree.size()), l(Distance_larger(Search_nearest)), max_k(k), distance_instance(d) 

  {
    if (search_nearest) 
      distance_to_root = d.min_distance_to_rectangle(q, tree.bounding_box());
    else 
      distance_to_root = d.max_distance_to_rectangle(q, tree.bounding_box());

    compute_neighbors_orthogonally(tree.root(), distance_to_root);

  }


  // Print statistics of the k_neighbor search process.
  std::ostream& statistics (std::ostream& s) 
  {
    s << "K_Neighbor search statistics:" << std::endl;
    s << "Number of internal nodes visited:" 
      << number_of_internal_nodes_visited << std::endl;
    s << "Number of leaf nodes visited:" 
      << number_of_leaf_nodes_visited << std::endl;
    s << "Number of items visited:" 
      << number_of_items_visited << std::endl;
    return s;
  }



private:

  void compute_neighbors_orthogonally(Node_handle N, FT rd)
  {
    typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
    typename SearchTraits::Cartesian_const_iterator_d query_object_it = construct_it(query_object);
    if (!(N->is_leaf())) 
    {
      number_of_internal_nodes_visited++;
      int new_cut_dim=N->cutting_dimension();
      FT old_off, new_rd;
      FT new_off = *(query_object_it + new_cut_dim) - N->cutting_value();
      if ( ((new_off <  FT(0.0)) && (search_nearest)) ||
           ((new_off >= FT(0.0)) && (!search_nearest)) ) 
      {
        compute_neighbors_orthogonally(N->lower(),rd);
        if (search_nearest) {
          old_off= *(query_object_it + new_cut_dim) - N->low_value();
          if (old_off>FT(0.0)) 
            old_off=FT(0.0);
        }
        else {	
          old_off= *(query_object_it + new_cut_dim) - N->high_value();
          if (old_off<FT(0.0)) 
            old_off=FT(0.0);
        }
        new_rd = distance_instance.new_distance(rd,old_off,new_off,new_cut_dim);
        if (branch(new_rd)) 
          compute_neighbors_orthogonally(N->upper(), new_rd);                               
      }
      else // compute new distance
      {
        compute_neighbors_orthogonally(N->upper(),rd); 
        if (search_nearest) {
          old_off= N->high_value() - *(query_object_it + new_cut_dim);
          if (old_off>FT(0.0)) 
            old_off=FT(0.0);
        }
        else  {       
          old_off= N->low_value() - *(query_object_it + new_cut_dim);
          if (old_off<FT(0.0)) 
            old_off=FT(0.0);
        }  
        new_rd = distance_instance. new_distance(rd,old_off,new_off,new_cut_dim);
        if (branch(new_rd)) 
          compute_neighbors_orthogonally(N->lower(), new_rd);       
      }
    }
    else
    {
      // n is a leaf
      number_of_leaf_nodes_visited++;
      if (N->size() > 0)
      {
        for (Point_d_iterator it=N->begin(); it != N->end(); it++) 
        {
          number_of_items_visited++;
          FT distance_to_query_object=
            distance_instance.transformed_distance(query_object,**it);
          insert(*it,distance_to_query_object);
        }
      }
    }
  }

}; // class 

} // namespace CGAL

#endif  // CGAL_ORTHOGONAL_K_NEIGHBOR_SEARCH_H
