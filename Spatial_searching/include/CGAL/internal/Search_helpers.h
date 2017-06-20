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
// Author(s)     : Clement Jamin (clement.jamin.pro@gmail.com)

#ifndef CGAL_INTERNAL_SEARCH_HELPERS_H
#define CGAL_INTERNAL_SEARCH_HELPERS_H

#include <CGAL/Has_member.h>

#include <boost/mpl/has_xxx.hpp>

namespace CGAL {
namespace internal {

// Helper struct to know at compile-time if there is a cache of the points
// stored in the tree
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

CGAL_GENERATE_MEMBER_DETECTOR(transformed_distance_from_coordinates);
CGAL_GENERATE_MEMBER_DETECTOR(interruptible_transformed_distance);
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(has_Enable_points_cache, Enable_points_cache, false)


template <typename Distance>
class Distance_helper
{
  typedef typename Distance::FT         FT;
  typedef typename Distance::Point_d    Point;
  typedef typename Distance::Query_item Query_item;

public:

  Distance_helper(Distance const& distance)
    : m_distance(distance)
  {}

  // If transformed_distance_from_coordinates does not exist in `Distance`
  template <bool has_transformed_distance_from_coordinates = has_transformed_distance_from_coordinates<Distance>::value>
  FT
  transformed_distance_from_coordinates(
    const typename Query_item& q,
    Point const& p,
    typename std::vector<FT>::const_iterator it_coord_begin,
    typename std::vector<FT>::const_iterator it_coord_end)
  {
    return m_distance.transformed_distance(q, p);
  }
  // ... or if it exists
  template <>
  FT
  transformed_distance_from_coordinates<true>(
    const typename Query_item& q,
    Point const& p,
    typename std::vector<FT>::const_iterator it_coord_begin,
    typename std::vector<FT>::const_iterator it_coord_end)
  {
    return m_distance.transformed_distance_from_coordinates(q, it_coord_begin, it_coord_end);
  }

  // *** Version with cache ***
  // If interruptible_transformed_distance does not exist in `Distance`
  template <bool has_interruptible_distance_computation = has_interruptible_transformed_distance<Distance>::value>
  FT
  interruptible_transformed_distance(
    const typename Query_item& q,
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
  interruptible_transformed_distance<true>(
    const typename Query_item& q,
    Point const& p,
    typename std::vector<FT>::const_iterator it_coord_begin,
    typename std::vector<FT>::const_iterator it_coord_end,
    FT stop_if_geq_to_this)
  {
    return m_distance.interruptible_transformed_distance(
      q, it_coord_begin, it_coord_end, stop_if_geq_to_this);
  }

  // *** Version without cache ***
  // If interruptible_transformed_distance does not exist in `Distance`
  template <bool has_interruptible_distance_computation = has_interruptible_transformed_distance<Distance>::value>
  FT
  interruptible_transformed_distance(
    const typename Query_item& q,
    Point const& p,
    FT)
  {
    return m_distance.transformed_distance(q, p);
  }
  // ... or if it exists
  template <>
  FT
  interruptible_transformed_distance<true>(
    const typename Query_item& q,
    Point const& p,
    FT stop_if_geq_to_this)
  {
    typename SearchTraits::Construct_cartesian_const_iterator_d construct_it = m_tree.traits().construct_cartesian_const_iterator_d_object();
    return m_distance.interruptible_transformed_distance(
      q, construct_it(p), construct_it(p, 0), stop_if_geq_to_this);
  }

private:
  Distance const & m_distance;

}; // Distance_helper


}} // namespace CGAL::internal

#endif  // CGAL_INTERNAL_SEARCH_HELPERSS_H
