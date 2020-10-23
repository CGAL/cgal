// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_POINT_SET_NEIGHBORHOOD_H
#define CGAL_CLASSIFICATION_POINT_SET_NEIGHBORHOOD_H

#include <CGAL/license/Classification.h>

#include <vector>

#include <CGAL/boost/iterator/counting_iterator.hpp>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/centroid.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/squared_distance_3.h>

#include <CGAL/Classification/Image.h>


namespace CGAL {

namespace Classification {

  /*!
    \ingroup PkgClassificationPointSet

    \brief Class that precomputes spatial searching structures for an
    input point set and gives access to the local neighborhood of a
    point as a set of indices.

    It allows the user to generate models of `NeighborQuery` based on
    a fixed range neighborhood or on a fixed K number of neighbors. In
    addition, the spatial searching structures can be computed on a
    simplified version of the point set to allow for neighbor queries
    at a higher scale.

    \tparam GeomTraits is a model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type is
    `RandomAccessIterator` and its value type is the key type of
    `PointMap`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `GeomTraits::Point_3`.
  */
template <typename GeomTraits, typename PointRange, typename PointMap>
class Point_set_neighborhood
{
  using FT = typename GeomTraits::FT;
  using Point = typename GeomTraits::Point_3;

  class My_point_property_map{
    const PointRange* input;
    PointMap point_map;

  public:
    using value_type = typename boost::property_traits<PointMap>::value_type;
    using reference = typename boost::property_traits<PointMap>::reference;
    using key_type = std::uint32_t;
    using category = typename boost::property_traits<PointMap>::category;
    My_point_property_map () { }
    My_point_property_map (const PointRange *input, PointMap point_map)
      : input (input), point_map (point_map) { }
    friend reference get (const My_point_property_map& ppmap, key_type i)
    { return get(ppmap.point_map, *(ppmap.input->begin()+std::size_t(i))); }
  };

  using Search_traits_base = Search_traits_3<GeomTraits>;
  using Search_traits =  Search_traits_adapter<std::uint32_t, My_point_property_map, Search_traits_base>;
  using Splitter = Sliding_midpoint<Search_traits>;
  using Distance = Distance_adapter<std::uint32_t, My_point_property_map, Euclidean_distance<Search_traits_base> >;
  using Tree = Kd_tree<Search_traits, Splitter, Tag_true, Tag_true>;
  using Sphere = Fuzzy_sphere<Search_traits>;
  using Knn = Orthogonal_k_neighbor_search<Search_traits, Distance, Splitter, Tree>;

  std::shared_ptr<Tree> m_tree;
  Distance m_distance;

public:

  /*!
    Functor that computes the neighborhood of an input point with a
    fixed number of neighbors.

    \cgalModels CGAL::Classification::NeighborQuery

    \sa Point_set_neighborhood
  */
  class K_neighbor_query
  {
  public:
    using value_type = typename Point_set_neighborhood::Point; ///<
  private:
    const Point_set_neighborhood& neighborhood;
    unsigned int k;
  public:
    /*!
      \brief constructs a K neighbor query object.
      \param neighborhood point set neighborhood object.
      \param k number of neighbors per query.
    */
    K_neighbor_query (const Point_set_neighborhood& neighborhood, unsigned int k)
      : neighborhood (neighborhood), k(k) { }

    /// \cond SKIP_IN_MANUAL
    template <typename OutputIterator>
    OutputIterator operator() (const value_type& query, OutputIterator output) const
    {
      neighborhood.k_neighbors (query, k, output);
      return output;
    }
    /// \endcond
  };

  /*!
    Functor that computes the neighborhood of an input point defined
    as the points lying in a sphere of fixed radius centered at the
    input point.

    \cgalModels CGAL::Classification::NeighborQuery

    \sa Point_set_neighborhood
  */
  class Sphere_neighbor_query
  {
  public:
    using value_type = typename Point_set_neighborhood::Point; ///<
  private:
    const Point_set_neighborhood& neighborhood;
    float radius;
  public:
    /*!
      \brief constructs a range neighbor query object.
      \param neighborhood point set neighborhood object.
      \param radius radius of the neighbor query sphere.
    */
    Sphere_neighbor_query (const Point_set_neighborhood& neighborhood, float radius)
      : neighborhood (neighborhood), radius(radius) { }

    /// \cond SKIP_IN_MANUAL
    template <typename OutputIterator>
    OutputIterator operator() (const value_type& query, OutputIterator output) const
    {
      neighborhood.sphere_neighbors (query, radius, output);
      return output;
    }
    /// \endcond
  };

  /// \cond SKIP_IN_MANUAL
  friend class K_neighbor_query;
  friend class Sphere_neighbor_query;

  Point_set_neighborhood () { }
  /// \endcond

  /// \name Constructors
  /// @{

  /*!
    \brief constructs a neighborhood object based on the input range.

    \tparam ConcurrencyTag enables sequential versus parallel
    algorithm. Possible values are `Sequential_tag`, `Parallel_tag`,
    and `Parallel_if_available_tag`. If no tag is provided,
    `Parallel_if_available_tag` is used.

    \param input point range.
    \param point_map property map to access the input points.
  */
  template <typename ConcurrencyTag>
  Point_set_neighborhood (const PointRange& input,
                          PointMap point_map,
                          const ConcurrencyTag&)
    : m_tree (nullptr)
  {
    init<ConcurrencyTag> (input, point_map);
  }

  /// \cond SKIP_IN_MANUAL
  Point_set_neighborhood (const PointRange& input, PointMap point_map)
    : m_tree (nullptr)
  {
    init<Parallel_if_available_tag> (input, point_map);
  }

  template <typename ConcurrencyTag>
  void init (const PointRange& input, PointMap point_map)
  {
    My_point_property_map pmap (&input, point_map);
    m_tree = std::make_shared<Tree>
      (boost::counting_iterator<std::uint32_t> (0),
       boost::counting_iterator<std::uint32_t> (std::uint32_t(input.size())),
       Splitter(),
       Search_traits (pmap));
    m_distance = Distance (pmap);
    m_tree->template build<ConcurrencyTag>();
  }
  /// \endcond

  /*!
    \brief constructs a simplified neighborhood object based on the input range.

    This method first computes a simplified version of the input point
    set by voxelization: a 3D grid is defined and for each subset
    present in one cell, only the point closest to the centroid of
    this subset is used.

    \tparam ConcurrencyTag enables sequential versus parallel
    algorithm. Possible values are `Sequential_tag`, `Parallel_tag`,
    and `Parallel_if_available_tag`. If no tag is provided,
    `Parallel_if_available_tag` is used.

    \param input input range.
    \param point_map property map to access the input points.
    \param voxel_size size of the cells of the 3D grid used for simplification.
  */
  template <typename ConcurrencyTag>
  Point_set_neighborhood (const PointRange& input,
                          PointMap point_map,
                          float voxel_size,
                          const ConcurrencyTag&)
    : m_tree (nullptr)
  {
    init<ConcurrencyTag> (input, point_map, voxel_size);
  }

  /// \cond SKIP_IN_MANUAL
  Point_set_neighborhood (const PointRange& input,
                          PointMap point_map,
                          float voxel_size)
  {
    init<Parallel_if_available_tag> (input, point_map, voxel_size);
  }

  template <typename ConcurrencyTag>
  void init (const PointRange& input, PointMap point_map, float voxel_size)
  {
    // First, simplify
    std::vector<std::uint32_t> indices;
    My_point_property_map pmap (&input, point_map);
    voxelize_point_set(input.size(), indices, pmap, voxel_size);

    m_tree = std::make_shared<Tree> (indices.begin(), indices.end(),
                                     Splitter(),
                                     Search_traits (pmap));
    m_distance = Distance (pmap);
    m_tree->template build<ConcurrencyTag>();
  }
  /// \endcond

  /// @}

  /// \name Queries
  /// @{

  /*!
    \brief returns a neighbor query object with fixed number of neighbors `k`.
  */
  K_neighbor_query k_neighbor_query (const unsigned int k) const
  {
    return K_neighbor_query (*this, k);
  }

  /*!
    \brief returns a neighbor query object with fixed radius `radius`.
  */
  Sphere_neighbor_query sphere_neighbor_query (const float radius) const
  {
    return Sphere_neighbor_query (*this, radius);
  }

  /// @}

private:

  template <typename OutputIterator>
  void sphere_neighbors (const Point& query, const FT radius_neighbors, OutputIterator output) const
  {
    CGAL_assertion (m_tree != nullptr);
    Sphere fs (query, radius_neighbors, 0, m_tree->traits());
    m_tree->search (output, fs);
  }

  template <typename OutputIterator>
  void k_neighbors (const Point& query, const unsigned int k, OutputIterator output) const
  {
    CGAL_assertion (m_tree != nullptr);
    Knn search (*m_tree, query, k, 0, true, m_distance);
    for (typename Knn::iterator it = search.begin(); it != search.end(); ++ it)
      *(output ++) = it->first;
  }

  template <typename Map>
  void voxelize_point_set (std::size_t nb_pts, std::vector<std::uint32_t>& indices, Map point_map,
                           float voxel_size)
  {
    std::map<Point, std::vector<std::uint32_t> > grid;

    for (std::uint32_t i = 0; i < nb_pts; ++ i)
    {
      const Point& p = get(point_map, i);
      Point ref (std::floor(p.x() / voxel_size),
                 std::floor(p.y() / voxel_size),
                 std::floor(p.z() / voxel_size));
      typename std::map<Point, std::vector<std::uint32_t> >::iterator it;
      boost::tie (it, boost::tuples::ignore)
        = grid.insert (std::make_pair (ref, std::vector<std::uint32_t>()));
      it->second.push_back (i);
    }

    for (typename std::map<Point, std::vector<std::uint32_t> >::iterator
           it = grid.begin(); it != grid.end(); ++ it)
    {
      const std::vector<std::uint32_t>& pts = it->second;
      Point centroid = CGAL::centroid (CGAL::make_transform_iterator_from_property_map
                                       (pts.begin(), point_map),
                                       CGAL::make_transform_iterator_from_property_map
                                       (pts.end(), point_map));
      std::uint32_t chosen = 0;
      float min_dist = (std::numeric_limits<float>::max)();
      for (std::size_t i = 0; i < pts.size(); ++ i)
      {
        float dist = float(CGAL::squared_distance(get(point_map, pts[i]), centroid));
        if (dist < min_dist)
        {
          min_dist = dist;
          chosen = pts[i];
        }
      }
      indices.push_back (chosen);
    }
  }
};


}

}


#endif // CGAL_CLASSIFICATION_POINT_SET_POINT_SET_NEIGHBORHOOD_H
