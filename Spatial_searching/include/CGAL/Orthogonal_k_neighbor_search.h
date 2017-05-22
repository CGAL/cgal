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
// 
//
// Author(s)     : Gael Guennebaud (gael.guennebaud@inria.fr), Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef CGAL_ORTHOGONAL_K_NEIGHBOR_SEARCH_H
#define CGAL_ORTHOGONAL_K_NEIGHBOR_SEARCH_H

#include <CGAL/internal/K_neighbor_search.h>
#include <CGAL/Has_member.h>

#include <boost/mpl/has_xxx.hpp>

namespace CGAL {

template <typename Tree, bool has_enable_points_cache>
struct Has_points_cache;

template <typename Tree>
struct Has_points_cache<Tree, true>
{
  static const bool value = Tree::Enable_points_cache::value;
};

template <typename Tree>
struct Has_points_cache<Tree, false>
{
  static const bool value = false;
};

template <class SearchTraits, 
          class Distance= typename internal::Spatial_searching_default_distance<SearchTraits>::type,
          class Splitter= Sliding_midpoint<SearchTraits> ,
          class Tree= Kd_tree<SearchTraits, Splitter, Tag_true> >
class Orthogonal_k_neighbor_search: public internal::K_neighbor_search<SearchTraits,Distance,Splitter,Tree>
{
  typedef internal::K_neighbor_search<SearchTraits,Distance,Splitter,Tree> Base;
  typedef typename Tree::Point_d Point;

public:
  typedef typename Base::FT FT;

private:
  typename SearchTraits::Cartesian_const_iterator_d query_object_it;
  
  std::vector<FT> dists;
  int m_dim;
  Tree const& m_tree;

  CGAL_GENERATE_MEMBER_DETECTOR(transformed_distance_from_coordinates);
  CGAL_GENERATE_MEMBER_DETECTOR(interruptable_transformed_distance);
  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(has_Enable_points_cache, Enable_points_cache, false)

  // If transformed_distance_from_coordinates does not exist in `Distance`
  template <bool has_transformed_distance_from_coordinates = has_transformed_distance_from_coordinates<Distance>::value>
  FT
  transformed_distance_from_coordinates(
    const typename Base::Query_item& q,
    Point const& p,
    typename std::vector<FT>::const_iterator it_coord_begin,
    typename std::vector<FT>::const_iterator it_coord_end)
  {
    return this->distance_instance.transformed_distance(q, p);
  }
  // ... or if it exists
  template <>
  FT
    transformed_distance_from_coordinates<true>(
    const typename Base::Query_item& q,
    Point const& p,
    typename std::vector<FT>::const_iterator it_coord_begin,
    typename std::vector<FT>::const_iterator it_coord_end)
  {
    return this->distance_instance.transformed_distance_from_coordinates(q, it_coord_begin, it_coord_end);
  }

  // *** Version with cache ***
  // If interruptable_transformed_distance does not exist in `Distance`
  template <bool has_interruptable_distance_computation = has_interruptable_transformed_distance<Distance>::value>
  FT
  interruptable_transformed_distance(
    const typename Base::Query_item& q,
    Point const& p,
    typename std::vector<FT>::const_iterator it_coord_begin, 
    typename std::vector<FT>::const_iterator it_coord_end,
    FT)
  {
    return transformed_distance_from_coordinates(q, p, it_coord_begin, it_coord_end);
  }
  // ... or if it exists
  template <>
  FT
  interruptable_transformed_distance<true>(
    const typename Base::Query_item& q,
    Point const& p,
    typename std::vector<FT>::const_iterator it_coord_begin,
    typename std::vector<FT>::const_iterator it_coord_end,
    FT stop_if_geq_to_this)
  {
    return this->distance_instance.interruptable_transformed_distance(
      q, it_coord_begin, it_coord_end, stop_if_geq_to_this);
  }

  // *** Version without cache ***
  // If interruptable_transformed_distance does not exist in `Distance`
  template <bool has_interruptable_distance_computation = has_interruptable_transformed_distance<Distance>::value>
  FT
  interruptable_transformed_distance(
    const typename Base::Query_item& q,
    Point const& p,
    FT)
  {
    return this->distance_instance.transformed_distance(q, p);
  }
  // ... or if it exists
  template <>
  FT
  interruptable_transformed_distance<true>(
    const typename Base::Query_item& q,
    Point const& p,
    FT stop_if_geq_to_this)
  {
    typename SearchTraits::Construct_cartesian_const_iterator_d construct_it = m_tree.traits().construct_cartesian_const_iterator_d_object();
    return this->distance_instance.interruptable_transformed_distance(
      q, construct_it(p), construct_it(p, 0), stop_if_geq_to_this);
  }

public:

  Orthogonal_k_neighbor_search(const Tree& tree, const typename Base::Query_item& q,  
                               unsigned int k=1, FT Eps=FT(0.0), bool Search_nearest=true, const Distance& d=Distance(),bool sorted=true)
    : Base(q,k,Eps,Search_nearest,d), m_tree(tree)
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

  template<bool use_cache>
  void search_in_leaf(typename Tree::Leaf_node_const_handle node);

  // With cache
  template<>
  void search_in_leaf<true>(typename Tree::Leaf_node_const_handle node)
  {
    typename Tree::iterator it_node_point = node->begin(), it_node_point_end = node->end();
    typename std::vector<FT>::const_iterator cache_point_begin = m_tree.cache_begin() + m_dim*(it_node_point - m_tree.begin());
    for (; !this->queue.full() && it_node_point != it_node_point_end; ++it_node_point)
    {
      this->number_of_items_visited++;
          
      FT distance_to_query_object =
        transformed_distance_from_coordinates(this->query_object, *it_node_point, cache_point_begin, cache_point_begin + m_dim);
      this->queue.insert(std::make_pair(&(*it_node_point), distance_to_query_object));

      cache_point_begin += m_dim;
    }
    FT worst_dist = this->queue.top().second;
    for (; it_node_point != it_node_point_end; ++it_node_point)
    {
      this->number_of_items_visited++;

      FT distance_to_query_object =
        interruptable_transformed_distance(this->query_object, *it_node_point, cache_point_begin, cache_point_begin + m_dim, worst_dist);

      if (distance_to_query_object < worst_dist)
      {
        this->queue.insert(std::make_pair(&(*it_node_point), distance_to_query_object));
        worst_dist = this->queue.top().second;
      }

      cache_point_begin += m_dim;
    }
  }

  // Without cache
  template<>
  void search_in_leaf<false>(typename Tree::Leaf_node_const_handle node)
  {
    typename Tree::iterator it_node_point = node->begin(), it_node_point_end = node->end();
    for (; !this->queue.full() && it_node_point != it_node_point_end; ++it_node_point)
    {
      this->number_of_items_visited++;

      FT distance_to_query_object =
        this->distance_instance.transformed_distance(this->query_object, *it_node_point);
      this->queue.insert(std::make_pair(&(*it_node_point), distance_to_query_object));
    }
    FT worst_dist = this->queue.top().second;
    for (; it_node_point != it_node_point_end; ++it_node_point)
    {
      this->number_of_items_visited++;

      FT distance_to_query_object =
        interruptable_transformed_distance(this->query_object, *it_node_point, worst_dist);

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
        search_in_leaf<Has_points_cache<Tree, has_Enable_points_cache<Tree>::type::value>::value>(node);
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
        for (typename Tree::iterator it_node_point=node->begin(); it_node_point != node->end(); it_node_point++) 
        {
          this->number_of_items_visited++;
          FT distance_to_query_object=
            this->distance_instance.transformed_distance(this->query_object,*it_node_point);
          this->queue.insert(std::make_pair(&(*it_node_point),distance_to_query_object));
        }
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
          new_off = 2*val < node->upper_low_value()+node->upper_high_value() ?
                    val - node->upper_high_value():
                    val - node->upper_low_value();
          bestChild = node->lower();
          otherChild = node->upper();
      }
      else // compute new distance
      {
          new_off = 2*val < node->lower_low_value()+node->lower_high_value() ?
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

#endif  // CGAL_ORTHOGONAL_K_NEIGHBOR_SEARCH_H
