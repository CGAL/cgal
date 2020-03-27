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
// Author(s)     : Clement Jamin (clement.jamin.pro@gmail.com)

#ifndef CGAL_INTERNAL_SEARCH_HELPERS_H
#define CGAL_INTERNAL_SEARCH_HELPERS_H

#include <CGAL/license/Spatial_searching.h>

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
  typedef typename Tree::Enable_points_cache type;
  static const bool value = type::value;
};

template <typename Tree>
struct Has_points_cache<Tree, false>
{
  typedef Tag_false type;
  static const bool value = false;
};

CGAL_GENERATE_MEMBER_DETECTOR(transformed_distance_from_coordinates);
CGAL_GENERATE_MEMBER_DETECTOR(interruptible_transformed_distance);
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(has_Enable_points_cache, Enable_points_cache, false)



// If transformed_distance_from_coordinates does not exist in `Distance`
template <typename Distance, typename SearchTraits, bool has_transformed_distance_from_coordinates>
class Transformed_distance_from_coordinates
{
  typedef typename Distance::FT         FT;
  typedef typename Distance::Point_d    Point;
  typedef typename Distance::Query_item Query_item;

public:
  Transformed_distance_from_coordinates(Distance const& distance)
    : m_distance(distance)
  {}

  FT operator() (
    const Query_item& q,
    Point const& p,
    typename std::vector<FT>::const_iterator /*it_coord_begin*/,
    typename std::vector<FT>::const_iterator /*it_coord_end*/) const
  {
    return m_distance.transformed_distance(q, p);
  }

private:
  Distance const& m_distance;
};
// ... or if it exists
template <typename Distance, typename SearchTraits>
class Transformed_distance_from_coordinates<Distance, SearchTraits, true>
{
  typedef typename Distance::FT         FT;
  typedef typename Distance::Point_d    Point;
  typedef typename Distance::Query_item Query_item;

public:
  Transformed_distance_from_coordinates(Distance const& distance)
    : m_distance(distance)
  {}

  FT operator() (
    const Query_item& q,
    Point const& /*p*/,
    typename std::vector<FT>::const_iterator it_coord_begin,
    typename std::vector<FT>::const_iterator it_coord_end) const
  {
    return m_distance.transformed_distance_from_coordinates(q, it_coord_begin, it_coord_end);
  }

private:
  Distance const& m_distance;
};

// If interruptible_transformed_distance does not exist in `Distance`
template <typename Distance, typename SearchTraits, bool has_interruptible_transformed_distance>
class Interruptible_transformed_distance
{
  typedef typename Distance::FT         FT;
  typedef typename Distance::Point_d    Point;
  typedef typename Distance::Query_item Query_item;

public:
  typedef Transformed_distance_from_coordinates<Distance, SearchTraits,
    has_transformed_distance_from_coordinates<Distance>::value> Tdfc;

  Interruptible_transformed_distance(
    SearchTraits const&, Distance const& distance, Tdfc const& tdfc)
    : m_distance(distance), m_ref_to_tdfc(tdfc)
  {}

  FT operator() (
    const Query_item& q,
    Point const& p,
    FT) const
  {
    return m_distance.transformed_distance(q, p);
  }

  FT operator() (
    const Query_item& q,
    Point const& p,
    typename std::vector<FT>::const_iterator it_coord_begin,
    typename std::vector<FT>::const_iterator it_coord_end,
    FT) const
  {
    return m_ref_to_tdfc(q, p, it_coord_begin, it_coord_end);
  }

private:
  Distance const& m_distance;
  Tdfc const& m_ref_to_tdfc;
};
// ... or if it exists
template <typename Distance, typename SearchTraits>
class Interruptible_transformed_distance<Distance, SearchTraits, true>
{
  typedef typename Distance::FT         FT;
  typedef typename Distance::Point_d    Point;
  typedef typename Distance::Query_item Query_item;

public:
  typedef Transformed_distance_from_coordinates<Distance, SearchTraits,
    has_transformed_distance_from_coordinates<Distance>::value> Tdfc;

  Interruptible_transformed_distance(
    SearchTraits const& traits, Distance const& distance, Tdfc const&)
    : m_traits(traits), m_distance(distance)
  {}

  FT operator() (
    const Query_item& q,
    Point const& p,
    FT stop_if_geq_to_this) const
  {
    typename SearchTraits::Construct_cartesian_const_iterator_d construct_it =
      m_traits.construct_cartesian_const_iterator_d_object();
    return m_distance.interruptible_transformed_distance(
      q, construct_it(p), construct_it(p, 0), stop_if_geq_to_this);
  }

  FT operator() (
    const Query_item& q,
    Point const& /*p*/,
    typename std::vector<FT>::const_iterator it_coord_begin,
    typename std::vector<FT>::const_iterator it_coord_end,
    FT stop_if_geq_to_this) const
  {
    return m_distance.interruptible_transformed_distance(
      q, it_coord_begin, it_coord_end, stop_if_geq_to_this);
  }

private:
  SearchTraits const& m_traits;
  Distance const& m_distance;
};



template <typename Distance, typename SearchTraits>
class Distance_helper
{
  typedef typename Distance::FT         FT;
  typedef typename Distance::Point_d    Point;
  typedef typename Distance::Query_item Query_item;

public:

  Distance_helper(Distance const& distance, SearchTraits const& traits)
    : m_distance(distance), m_tdfc(m_distance), m_itd(traits, distance, m_tdfc)
  {}

  FT
  transformed_distance_from_coordinates(
    const Query_item& q,
    Point const& p,
    typename std::vector<FT>::const_iterator it_coord_begin,
    typename std::vector<FT>::const_iterator it_coord_end)
  {
    return m_tdfc(q, p, it_coord_begin, it_coord_end);
  }

  FT
  interruptible_transformed_distance(
    const Query_item& q,
    Point const& p,
    FT stop_if_geq_to_this)
  {
    return m_itd(q, p, stop_if_geq_to_this);
  }

  FT
  interruptible_transformed_distance(
    const Query_item& q,
    Point const& p,
    typename std::vector<FT>::const_iterator it_coord_begin,
    typename std::vector<FT>::const_iterator it_coord_end,
    FT stop_if_geq_to_this)
  {
    return m_itd(q, p, it_coord_begin, it_coord_end, stop_if_geq_to_this);
  }

private:
  Distance     const& m_distance;
  Transformed_distance_from_coordinates<Distance, SearchTraits,
    has_transformed_distance_from_coordinates<Distance>::value> m_tdfc;
  Interruptible_transformed_distance<Distance, SearchTraits,
    has_interruptible_transformed_distance<Distance>::value> m_itd;
}; // Distance_helper


}} // namespace CGAL::internal

#endif  // CGAL_INTERNAL_SEARCH_HELPERSS_H
