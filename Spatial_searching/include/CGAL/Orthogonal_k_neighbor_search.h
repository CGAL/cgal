// Copyright (c) 2002,2011 Utrecht University (The Netherlands).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Gael Guennebaud (gael.guennebaud@inria.fr), Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef CGAL_ORTHOGONAL_K_NEIGHBOR_SEARCH_H
#define CGAL_ORTHOGONAL_K_NEIGHBOR_SEARCH_H

#include <CGAL/license/Spatial_searching.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/internal/K_neighbor_search.h>

namespace CGAL {

template <class SearchTraits, 
          class Distance= typename internal::Spatial_searching_default_distance<SearchTraits>::type,
          class Splitter= Sliding_midpoint<SearchTraits> ,
          class Tree= Kd_tree<SearchTraits, Splitter, Tag_true> >
class Orthogonal_k_neighbor_search: public internal::K_neighbor_search<SearchTraits,Distance,Splitter,Tree> {
  typedef  internal::K_neighbor_search<SearchTraits,Distance,Splitter,Tree> Base;

  typename SearchTraits::Cartesian_const_iterator_d query_object_it;
  
  std::vector<typename Base::FT> dists;
public:
  typedef typename Base::FT FT;


  Orthogonal_k_neighbor_search(const Tree& tree, const typename Base::Query_item& q,  
                               unsigned int k=1, FT Eps=FT(0.0), bool Search_nearest=true, const Distance& d=Distance(),bool sorted=true)
    : Base(q,k,Eps,Search_nearest,d) 
  {
    if (tree.empty()) return;

    typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=tree.traits().construct_cartesian_const_iterator_d_object();
    query_object_it = construct_it(this->query_object);

    int dim = static_cast<int>(std::distance(query_object_it, construct_it(this->query_object,0)));

    dists.resize(dim);
    for(int i=0;i<dim;i++)
        dists[i]=0;

    FT distance_to_root;
    if (this->search_nearest){ 
      distance_to_root = this->distance_instance.min_distance_to_rectangle(q, tree.bounding_box(),dists);
      compute_nearest_neighbors_orthogonally(tree.root(), distance_to_root);
    }
    else {
      distance_to_root = this->distance_instance.max_distance_to_rectangle(q, tree.bounding_box(),dists);
      compute_furthest_neighbors_orthogonally(tree.root(), distance_to_root);
    }

    
      
          
      
    if (sorted) this->queue.sort();
  }
private:

  void compute_nearest_neighbors_orthogonally(typename Base::Node_const_handle N, FT rd)
  {
    if (!(N->is_leaf())) 
    {
      typename Tree::Internal_node_const_handle node =
        static_cast<typename Tree::Internal_node_const_handle>(N);
      this->number_of_internal_nodes_visited++;
      int new_cut_dim=node->cutting_dimension();
      typename Base::Node_const_handle bestChild, otherChild;
      FT new_off;
      FT val = *(query_object_it + new_cut_dim);
      FT diff1 = val - node->upper_low_value();
      FT diff2 = val - node->lower_high_value();
      if ( (diff1 + diff2 <  FT(0.0)) ) 
      {
          new_off = diff1;
          bestChild = node->lower();
          otherChild = node->upper();                      
      }
      else // compute new distance
      {
          new_off= diff2;
          bestChild = node->upper();
          otherChild = node->lower();
      }
      compute_nearest_neighbors_orthogonally(bestChild, rd);
      FT dst=dists[new_cut_dim];
      FT new_rd = this->distance_instance.new_distance(rd,dst,new_off,new_cut_dim);
      dists[new_cut_dim]=new_off;
        if (this->branch_nearest(new_rd)) 
        {
          compute_nearest_neighbors_orthogonally(otherChild, new_rd);
        }
      dists[new_cut_dim]=dst;
    }
    else
    {
      // n is a leaf
      typename Tree::Leaf_node_const_handle node =
        static_cast<typename Tree::Leaf_node_const_handle>(N);
      this->number_of_leaf_nodes_visited++;
      bool full = this->queue.full();
      FT worst_dist = this->queue.top().second;
      if (node->size() > 0)
      {
        for (typename Tree::iterator it=node->begin(); it != node->end(); it++) 
        {
          this->number_of_items_visited++;
          FT distance_to_query_object=
            this->distance_instance.transformed_distance(this->query_object,*it);
          
          if(!full || distance_to_query_object < worst_dist)
            this->queue.insert(std::make_pair(&(*it),distance_to_query_object));
        }
      }
    }
  }    

   void compute_furthest_neighbors_orthogonally(typename Base::Node_const_handle N, FT rd)
  {
    if (!(N->is_leaf())) 
    {
      typename Tree::Internal_node_const_handle node = 
        static_cast<typename Tree::Internal_node_const_handle>(N);
      this->number_of_internal_nodes_visited++;
      int new_cut_dim=node->cutting_dimension();
      typename Base::Node_const_handle bestChild, otherChild;
      FT new_off;
      FT val = *(query_object_it + new_cut_dim);
      FT diff1 = val - node->lower_high_value();
      FT diff2 = val - node->upper_low_value();
      if ( (diff1 + diff2 >= FT(0.0)) ) 
      {
          new_off = node->upper_low_value()+node->upper_high_value() > val*2?
                    val - node->upper_high_value():
                    val - node->upper_low_value();
          bestChild = node->lower();
          otherChild = node->upper();
      }
      else // compute new distance
      {
          new_off = node->lower_low_value()+node->lower_high_value() > val*2 ?
                    val - node->lower_high_value():
                    val - node->lower_low_value();
          bestChild = node->upper();
          otherChild = node->lower();
      }
      compute_furthest_neighbors_orthogonally(bestChild,rd);
      FT dst=dists[new_cut_dim];
      FT new_rd = this->distance_instance.new_distance(rd,dst,new_off,new_cut_dim);
      dists[new_cut_dim]=new_off;
        if (this->branch_furthest(new_rd)) 
          compute_furthest_neighbors_orthogonally(otherChild, new_rd);
      dists[new_cut_dim]=dst;
    }
    else
    {
      // n is a leaf
      typename Tree::Leaf_node_const_handle node = 
        static_cast<typename Tree::Leaf_node_const_handle>(N);
      this->number_of_leaf_nodes_visited++;
      if (node->size() > 0)
      {
        for (typename Tree::iterator it=node->begin(); it != node->end(); it++) 
        {
          this->number_of_items_visited++;
          FT distance_to_query_object=
            this->distance_instance.transformed_distance(this->query_object,*it);
          this->queue.insert(std::make_pair(&(*it),distance_to_query_object));
        }
      }
    }
  }    
  
}; // class 

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif  // CGAL_ORTHOGONAL_K_NEIGHBOR_SEARCH_H
