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
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>),
//                 Clement Jamin (clement.jamin.pro@gmail.com)

#ifndef  CGAL_K_NEIGHBOR_SEARCH_H
#define  CGAL_K_NEIGHBOR_SEARCH_H

#include <CGAL/license/Spatial_searching.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/internal/K_neighbor_search.h>
#include <CGAL/internal/Get_dimension_tag.h>
#include <CGAL/internal/Search_helpers.h>

#include <iterator> // for std::distance


namespace CGAL {
        template <class SearchTraits, class Distance,class Splitter,class Tree>
class K_neighbor_search;

template <class SearchTraits,
          class Distance= typename internal::Spatial_searching_default_distance<SearchTraits>::type,
          class Splitter= Sliding_midpoint<SearchTraits> ,
          class Tree= Kd_tree<SearchTraits, Splitter, Tag_true, Tag_false> >
class K_neighbor_search: public internal::K_neighbor_search<SearchTraits,Distance,Splitter,Tree>
{
  typedef  internal::K_neighbor_search<SearchTraits,Distance,Splitter,Tree> Base;
  typedef typename Tree::Point_d Point;

public:
  typedef typename Base::FT FT;
  typedef typename internal::Get_dimension_tag<SearchTraits>::Dimension D;

  K_neighbor_search(const Tree& tree, const typename Base::Query_item& q,
    unsigned int k=1, FT Eps=FT(0.0), bool Search_nearest=true, const Distance& d=Distance(),bool sorted=true)
  : Base(q,k,Eps,Search_nearest,d),
    m_distance_helper(this->distance_instance, tree.traits()),
    m_tree(tree)
  {
    if (tree.empty()) return;

    compute_neighbors_general(tree.root(),tree.bounding_box());
    if (sorted) this->queue.sort();
  };

private:
  typedef typename Base::Node_const_handle Node_const_handle;
  using Base::branch;

  internal::Distance_helper<Distance, SearchTraits> m_distance_helper;
  Tree const& m_tree;

  void
  compute_neighbors_general(typename Base::Node_const_handle N, const Kd_tree_rectangle<FT,D>& r)
  {
    if (!(N->is_leaf())) {
       typename Tree::Internal_node_const_handle node =
        static_cast<typename Tree::Internal_node_const_handle>(N);
      this->number_of_internal_nodes_visited++;
      int new_cut_dim=node->cutting_dimension();
      FT  new_cut_val=node->cutting_value();

      Kd_tree_rectangle<FT,D> r_lower(r);

      // modifies also r_lower to lower half
      Kd_tree_rectangle<FT,D> r_upper(r_lower);
      r_lower.split(r_upper, new_cut_dim, new_cut_val);

      FT distance_to_lower_half;
      FT distance_to_upper_half;

      if (this->search_nearest) {

        distance_to_lower_half =
          this->distance_instance. min_distance_to_rectangle(this->query_object,
                                                       r_lower);

        distance_to_upper_half =
          this->distance_instance.min_distance_to_rectangle(this->query_object,
                                                      r_upper);


      }
      else
        {

          distance_to_lower_half =
            this->distance_instance.max_distance_to_rectangle(this->query_object,
                                                        r_lower);

          distance_to_upper_half =
            this->distance_instance.max_distance_to_rectangle(this->query_object,
                                                        r_upper);

        }

      if ( (( this->search_nearest) &&
            (distance_to_lower_half < distance_to_upper_half))
           ||
           ((!this->search_nearest) &&
            (distance_to_lower_half >=
             distance_to_upper_half))  )
        {
          compute_neighbors_general(node->lower(), r_lower);
          if (branch(distance_to_upper_half))
            compute_neighbors_general (node->upper(), r_upper);
        }
      else
        {        compute_neighbors_general(node->upper(), r_upper);
        if (branch(distance_to_lower_half))
          compute_neighbors_general (node->lower(),
                                     r_lower);
        }

    }
    else
    {
      // n is a leaf
      typename Tree::Leaf_node_const_handle node =
        static_cast<typename Tree::Leaf_node_const_handle>(N);
      this->number_of_leaf_nodes_visited++;
      if (node->size() > 0)
      {
        typename internal::Has_points_cache<Tree, internal::has_Enable_points_cache<Tree>::type::value>::type dummy;
        if (this->search_nearest)
          search_nearest_in_leaf(node, dummy);
        else
          search_furthest_in_leaf(node, dummy);
      }
    }
  }

  // With cache
  void search_nearest_in_leaf(typename Tree::Leaf_node_const_handle node, Tag_true)
  {
    typename Tree::iterator it_node_point = node->begin(), it_node_point_end = node->end();
    int dim = m_tree.dim();
    typename std::vector<FT>::const_iterator cache_point_begin = m_tree.cache_begin() + dim*(it_node_point - m_tree.begin());
    // As long as the queue is not full, the node is just inserted
    for (; !this->queue.full() && it_node_point != it_node_point_end; ++it_node_point)
    {
      this->number_of_items_visited++;

      FT distance_to_query_object =
        m_distance_helper.transformed_distance_from_coordinates(
          this->query_object, *it_node_point, cache_point_begin, cache_point_begin + dim);
      this->queue.insert(std::make_pair(&(*it_node_point), distance_to_query_object));

      cache_point_begin += dim;
    }
    // Now that the queue is full, we can gain time by keeping track
    // of the current worst distance to interrupt the distance computation earlier
    FT worst_dist = this->queue.top().second;
    for (; it_node_point != it_node_point_end; ++it_node_point)
    {
      this->number_of_items_visited++;

      FT distance_to_query_object =
        m_distance_helper.interruptible_transformed_distance(
          this->query_object, *it_node_point, cache_point_begin, cache_point_begin + dim, worst_dist);

      if (distance_to_query_object < worst_dist)
      {
        this->queue.insert(std::make_pair(&(*it_node_point), distance_to_query_object));
        worst_dist = this->queue.top().second;
      }

      cache_point_begin += dim;
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


  // With cache
  void search_furthest_in_leaf(typename Tree::Leaf_node_const_handle node, Tag_true)
  {
    typename Tree::iterator it_node_point = node->begin(), it_node_point_end = node->end();
    int dim = m_tree.dim();
    typename std::vector<FT>::const_iterator cache_point_begin = m_tree.cache_begin() + dim*(it_node_point - m_tree.begin());
    // In furthest search mode, the interruptible distance cannot be used to optimize
    for (; it_node_point != it_node_point_end; ++it_node_point)
    {
      this->number_of_items_visited++;

      FT distance_to_query_object =
        m_distance_helper.transformed_distance_from_coordinates(
          this->query_object, *it_node_point, cache_point_begin, cache_point_begin + dim);
      this->queue.insert(std::make_pair(&(*it_node_point), distance_to_query_object));

      cache_point_begin += dim;
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

}; // class

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif  // CGAL_K_NEIGHBOR_SEARCH_H
