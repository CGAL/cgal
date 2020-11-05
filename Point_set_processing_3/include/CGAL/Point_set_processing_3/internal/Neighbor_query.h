// Copyright (c) 2019 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_PSP_INTERNAL_NEIGHBOR_QUERY_H
#define CGAL_PSP_INTERNAL_NEIGHBOR_QUERY_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/point_set_processing_assertions.h>

#include <CGAL/iterator.h>

#include <boost/function_output_iterator.hpp>

namespace CGAL {
namespace Point_set_processing_3 {
namespace internal {

struct Maximum_points_reached_exception : public std::exception { };

template <typename Kernel_, typename PointRangeRef, typename PointMap>
class Neighbor_query
{
public:

  typedef Kernel_ Kernel;
  typedef PointRangeRef Point_range;
  typedef PointMap Point_map;

  typedef typename Kernel::FT FT;
  typedef typename boost::property_traits<PointMap>::value_type Point;

  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Point_3 Point_3;

  typedef std::is_same<Point, Point_2> Is_2d;

  typedef typename Range_iterator_type<PointRangeRef>::type input_iterator;
  typedef typename input_iterator::value_type value_type;

  typedef CGAL::Prevent_deref<input_iterator> iterator;

  struct Deref_point_map
  {
    typedef input_iterator key_type;
    typedef typename boost::property_traits<PointMap>::reference reference;
    typedef typename boost::property_traits<PointMap>::value_type value_type;
    typedef typename boost::property_traits<PointMap>::category category;

    PointMap point_map;

    Deref_point_map () { }
    Deref_point_map (PointMap point_map) : point_map(point_map) { }

    friend reference get (const Deref_point_map& map, key_type it)
    {
      return get(map.point_map, *it);
    }
  };

  typedef typename std::conditional<Is_2d::value,
                                    CGAL::Search_traits_2<Kernel>,
                                    CGAL::Search_traits_3<Kernel> >::type Tree_traits_base;

  typedef CGAL::Search_traits_adapter<input_iterator, Deref_point_map, Tree_traits_base> Tree_traits;
  typedef CGAL::Sliding_midpoint<Tree_traits> Splitter;
  typedef CGAL::Distance_adapter<input_iterator, Deref_point_map, CGAL::Euclidean_distance<Tree_traits_base> > Distance;
  typedef CGAL::Kd_tree<Tree_traits, Splitter, CGAL::Tag_true, CGAL::Tag_true> Tree;
  typedef CGAL::Fuzzy_sphere<Tree_traits> Sphere;

  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits, Distance, Splitter, Tree> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

private:

  PointRangeRef m_points;
  PointMap m_point_map;
  Deref_point_map m_deref_map;
  Tree_traits m_traits;
  Tree m_tree;
  Distance m_distance;

  // Forbid copy
  Neighbor_query (const Neighbor_query&) { }

public:

  Neighbor_query (PointRangeRef points, PointMap point_map)
    : m_points (points)
    , m_point_map (point_map)
    , m_deref_map (point_map)
    , m_traits (m_deref_map)
    , m_tree (iterator(m_points.begin()), iterator(m_points.end()), Splitter(), m_traits)
    , m_distance (m_deref_map)
  {
    m_tree.build();
  }

  PointMap point_map() const { return m_point_map; }

  template <typename OutputIterator>
  void get_iterators (const Point& query, unsigned int k, FT neighbor_radius,
                      OutputIterator output, unsigned int fallback_k_if_sphere_empty = 3) const
  {
    if (neighbor_radius != FT(0))
    {
      Sphere fs (query, neighbor_radius, 0, m_traits);

      // if k=0, no limit on the number of neighbors returned
      if (k == 0)
        k = (std::numeric_limits<unsigned int>::max)();

      unsigned int nb = 0;

      try
      {
        std::function<void(const input_iterator&)> output_iterator_with_limit
          = [&](const input_iterator& it)
          {
            *(output ++) = it;
            if (++ nb == k)
              throw Maximum_points_reached_exception();
          };

        auto function_output_iterator
          = boost::make_function_output_iterator (output_iterator_with_limit);

        m_tree.search (function_output_iterator, fs);
      }
      catch (const Maximum_points_reached_exception&)
      { }

      // Fallback, if not enough  points are return, search for the knn
      // first points
      if (nb < fallback_k_if_sphere_empty)
        k = fallback_k_if_sphere_empty;
      else
        k = 0;
    }

    if (k != 0)
    {
      // Gather set of (k+1) neighboring points.
      // Perform k+1 queries (as in point set, the query point is
      // output first). Search may be aborted if k is greater
      // than number of input points.

      Neighbor_search search (m_tree, query, k+1, 0, true, m_distance);
      Search_iterator search_iterator = search.begin();
      unsigned int i;
      for (i = 0; i < (k+1); ++ i)
      {
        if(search_iterator == search.end())
          break; // premature ending
        *(output ++) = search_iterator->first;
        search_iterator++;
      }
    }
  }

  template <typename OutputIterator>
  void get_points (const Point& query, unsigned int k, FT neighbor_radius,
                   OutputIterator output, unsigned int fallback_k_if_sphere_empty = 3) const
  {
    return get_iterators(query, k, neighbor_radius,
                         boost::make_function_output_iterator
                         ([&](const input_iterator& it)
                          {
                            *(output ++) = get (m_point_map, *it);
                          }), fallback_k_if_sphere_empty);
  }
};

} } } // namespace CGAL::Point_set_processing_3::internal

#endif // CGAL_PSP_INTERNAL_NEIGHBOR_QUERY_H
