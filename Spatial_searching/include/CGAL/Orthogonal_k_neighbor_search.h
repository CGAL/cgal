// Copyright (c) 2002,2011 Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Gael Guennebaud (gael.guennebaud@inria.fr),
//                 Hans Tangelder (<hanst@cs.uu.nl>),
//                 Clement Jamin (clement.jamin.pro@gmail.com)

#ifndef CGAL_ORTHOGONAL_K_NEIGHBOR_SEARCH_H
#define CGAL_ORTHOGONAL_K_NEIGHBOR_SEARCH_H

#include <CGAL/license/Spatial_searching.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/internal/K_neighbor_search.h>
#include <CGAL/internal/Search_helpers.h>

#include <iterator> // for std::distance

namespace CGAL {

template <class SearchTraits,
          class Distance= typename internal::Spatial_searching_default_distance<SearchTraits>::type,
          class Splitter= Sliding_midpoint<SearchTraits> ,
          class Tree= Kd_tree<SearchTraits, Splitter, Tag_true, Tag_false> >
class Orthogonal_k_neighbor_search: public internal::K_neighbor_search<SearchTraits,Distance,Splitter,Tree>
{
  typedef internal::K_neighbor_search<SearchTraits,Distance,Splitter,Tree> Base;
  typedef typename Tree::Point_d Point;

public:
  typedef typename Base::FT FT;

private:
  typename SearchTraits::Cartesian_const_iterator_d query_object_it;

  internal::Distance_helper<Distance, SearchTraits> m_distance_helper;
  std::vector<FT> dists;
  int m_dim;
  Tree const& m_tree;

public:

  Orthogonal_k_neighbor_search(const Tree& tree, const typename Base::Query_item& q,
                               unsigned int k=1, FT Eps=FT(0.0), bool Search_nearest=true, const Distance& d=Distance(),bool sorted=true)
  : Base(q,k,Eps,Search_nearest,d),
    m_distance_helper(this->distance_instance, tree.traits()),
    m_tree(tree)
  {
    if (tree.empty()) return;

    typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=tree.traits().construct_cartesian_const_iterator_d_object();
    query_object_it = construct_it(this->query_object);

    m_dim = static_cast<int>(std::distance(query_object_it, construct_it(this->query_object,0)));

    dists.resize(m_dim);
    for(int i=0;i<m_dim;i++)
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

  // With cache
  void search_nearest_in_leaf(typename Tree::Leaf_node_const_handle node, Tag_true)
  {
    typename Tree::iterator it_node_point = node->begin(), it_node_point_end = node->end();
    typename std::vector<FT>::const_iterator cache_point_begin = m_tree.cache_begin() + m_dim*(it_node_point - m_tree.begin());
    // As long as the queue is not full, the node is just inserted
    for (; !this->queue.full() && it_node_point != it_node_point_end; ++it_node_point)
    {
      this->number_of_items_visited++;

      FT distance_to_query_object =
        m_distance_helper.transformed_distance_from_coordinates(
          this->query_object, *it_node_point, cache_point_begin, cache_point_begin + m_dim);
      this->queue.insert(std::make_pair(&(*it_node_point), distance_to_query_object));

      cache_point_begin += m_dim;
    }
    // Now that the queue is full, we can gain time by keeping track
    // of the current worst distance to interrupt the distance computation earlier
    FT worst_dist = this->queue.top().second;
    for (; it_node_point != it_node_point_end; ++it_node_point)
    {
      this->number_of_items_visited++;

      FT distance_to_query_object =
        m_distance_helper.interruptible_transformed_distance(
          this->query_object, *it_node_point, cache_point_begin, cache_point_begin + m_dim, worst_dist);

      if (distance_to_query_object < worst_dist)
      {
        this->queue.insert(std::make_pair(&(*it_node_point), distance_to_query_object));
        worst_dist = this->queue.top().second;
      }

      cache_point_begin += m_dim;
    }
  }

  // Without cache
  void search_nearest_in_leaf(typename Tree::Leaf_node_const_handle node, Tag_false)
  {
    typename Tree::iterator it_node_point = node->begin(), it_node_point_end = node->end();
    // As long as the queue is not full, the node is just inserted
    for (; !this->queue.full() && it_node_point != it_node_point_end; ++it_node_point)
    {
      this->number_of_items_visited++;

      FT distance_to_query_object =
        this->distance_instance.transformed_distance(this->query_object, *it_node_point);
      this->queue.insert(std::make_pair(&(*it_node_point), distance_to_query_object));
    }
    // Now that the queue is full, we can gain time by keeping track
    // of the current worst distance to interrupt the distance computation earlier
    FT worst_dist = this->queue.top().second;
    for (; it_node_point != it_node_point_end; ++it_node_point)
    {
      this->number_of_items_visited++;

      FT distance_to_query_object =
        m_distance_helper.interruptible_transformed_distance(
          this->query_object, *it_node_point, worst_dist);

      if (distance_to_query_object < worst_dist)
      {
        this->queue.insert(std::make_pair(&(*it_node_point), distance_to_query_object));
        worst_dist = this->queue.top().second;
      }
    }
  }

  void compute_nearest_neighbors_orthogonally(typename Base::Node_const_handle N, FT rd)
  {
    if (N->is_leaf())
    {
      // n is a leaf
      typename Tree::Leaf_node_const_handle node =
        static_cast<typename Tree::Leaf_node_const_handle>(N);
      this->number_of_leaf_nodes_visited++;
      if (node->size() > 0)
      {
        typename internal::Has_points_cache<Tree, internal::has_Enable_points_cache<Tree>::type::value>::type dummy;
        search_nearest_in_leaf(node, dummy);
      }
    }
    else
    {
      typename Tree::Internal_node_const_handle node =
        static_cast<typename Tree::Internal_node_const_handle>(N);
      this->number_of_internal_nodes_visited++;
      int new_cut_dim = node->cutting_dimension();
      typename Base::Node_const_handle bestChild, otherChild;
      FT new_off;
      FT val = *(query_object_it + new_cut_dim);
      FT diff1 = val - node->upper_low_value();
      FT diff2 = val - node->lower_high_value();
      if ((diff1 + diff2 <  FT(0.0)))
      {
        new_off = diff1;
        bestChild = node->lower();
        otherChild = node->upper();
      }
      else // compute new distance
      {
        new_off = diff2;
        bestChild = node->upper();
        otherChild = node->lower();
      }
      compute_nearest_neighbors_orthogonally(bestChild, rd);
      FT dst = dists[new_cut_dim];
      FT new_rd = this->distance_instance.new_distance(rd, dst, new_off, new_cut_dim);
      dists[new_cut_dim] = new_off;
      if (this->branch_nearest(new_rd))
      {
        compute_nearest_neighbors_orthogonally(otherChild, new_rd);
      }
      dists[new_cut_dim] = dst;
    }
  }

  // With cache
  void search_furthest_in_leaf(typename Tree::Leaf_node_const_handle node, Tag_true)
  {
    typename Tree::iterator it_node_point = node->begin(), it_node_point_end = node->end();
    typename std::vector<FT>::const_iterator cache_point_begin = m_tree.cache_begin() + m_dim*(it_node_point - m_tree.begin());
    // In furthest search mode, the interruptible distance cannot be used to optimize
    for (; it_node_point != it_node_point_end; ++it_node_point)
    {
      this->number_of_items_visited++;

      FT distance_to_query_object =
        m_distance_helper.transformed_distance_from_coordinates(
          this->query_object, *it_node_point, cache_point_begin, cache_point_begin + m_dim);
      this->queue.insert(std::make_pair(&(*it_node_point), distance_to_query_object));

      cache_point_begin += m_dim;
    }
  }

  // Without cache
  void search_furthest_in_leaf(typename Tree::Leaf_node_const_handle node, Tag_false)
  {
    typename Tree::iterator it_node_point = node->begin(), it_node_point_end = node->end();
    // In furthest search mode, the interruptible distance cannot be used to optimize
    for (; it_node_point != it_node_point_end; ++it_node_point)
    {
      this->number_of_items_visited++;

      FT distance_to_query_object =
        this->distance_instance.transformed_distance(this->query_object, *it_node_point);
      this->queue.insert(std::make_pair(&(*it_node_point), distance_to_query_object));
    }
  }

  void compute_furthest_neighbors_orthogonally(typename Base::Node_const_handle N, FT rd)
  {
    if (N->is_leaf())
    {
      // n is a leaf
      typename Tree::Leaf_node_const_handle node =
        static_cast<typename Tree::Leaf_node_const_handle>(N);
      this->number_of_leaf_nodes_visited++;
      if (node->size() > 0)
      {
        typename internal::Has_points_cache<Tree, internal::has_Enable_points_cache<Tree>::type::value>::type dummy;
        search_furthest_in_leaf(node, dummy);
      }
    }
    else
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
  }

}; // class

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif  // CGAL_ORTHOGONAL_K_NEIGHBOR_SEARCH_H
