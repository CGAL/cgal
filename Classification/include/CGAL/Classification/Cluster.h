// Copyright (c) 2018 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_CLUSTER_H
#define CGAL_CLASSIFICATION_CLUSTER_H

#include <CGAL/license/Classification.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/property_map.h>

#include <boost/iterator/transform_iterator.hpp>

namespace CGAL {

namespace Classification {


  /*!
    \ingroup PkgClassificationCluster

    \brief Class that represent a cluster of items to be classified as
    a single atomic object.

    A cluster is a set of indices of items inside an input range with
    random access.

    \tparam ItemRange model of `ConstRange`. Its iterator type is
    `RandomAccessIterator`. Its value type depends on the data that is
    classified (for example, `CGAL::Point_3` or `CGAL::Triangle_3`).

    \tparam ItemMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `ItemRange` and value type
    is the type of item to classify (for example, `CGAL::Point_3`).
  */
template <typename ItemRange, typename ItemMap>
class Cluster
{
public:

  using Item = typename boost::property_traits<ItemMap>::value_type;

  /// \cond SKIP_IN_MANUAL
  struct Neighbor_query
  {
    template <typename OutputIterator>
    OutputIterator operator() (const Cluster& cluster, OutputIterator output) const
    {
      return std::copy (cluster.neighbors->begin(), cluster.neighbors->end(), output);
    }
  };
  std::shared_ptr<std::vector<std::size_t> > neighbors;
  /// \endcond

  /// \cond SKIP_IN_MANUAL
  class Point_idx_to_point_unary_function
  {
  public:
    using argument_type = std::size_t;
    using result_type = typename boost::property_traits<ItemMap>::reference;
    using category = boost::readable_property_map_tag;

    const ItemRange* m_range;
    ItemMap m_item_map;

    Point_idx_to_point_unary_function (const ItemRange* range, ItemMap item_map)
      : m_range (range), m_item_map (item_map)
    { }

    result_type operator() (const argument_type& arg) const
    {
      return get (m_item_map, *(m_range->begin() + arg));
    }
  };
  /// \endcond

private:
  const ItemRange* m_range;
  ItemMap m_item_map;

  std::shared_ptr<std::vector<std::size_t> > m_inliers;
  mutable CGAL::Bbox_3 m_bounding_box;
  int m_training;
  int m_label;

public:

  /// \name Constructor
  /// @{

  /*!
    \brief Constructs an empty cluster of items.

    Items in the clusters will be subsets of `range`.

    \param range input range.
    \param item_map property map to access the input items.
  */
  Cluster (const ItemRange& range, ItemMap item_map)
    : neighbors (std::make_shared<std::vector<std::size_t> >())
    , m_range (&range), m_item_map (item_map)
    , m_inliers (std::make_shared<std::vector<std::size_t> >())
    , m_training(-1), m_label(-1)
  { }

  /// @}

  /// \name Modifications
  /// @{

  /*!
    \brief Clears the cluster.
  */
  void clear () { m_inliers->clear(); }

  /*!
    \brief Inserts element of index `idx` in the cluster.
  */
  void insert (std::size_t idx) { m_inliers->push_back (idx); }

  /// @}

  /// \name Access
  /// @{

  /*!
    \brief Returns the number of items in the cluster.
  */
  std::size_t size() const { return m_inliers->size(); }

  /*!
    \brief Returns the index (in the input range) of the i^{th} element of the cluster.
  */
  std::size_t index (std::size_t i) const { return (*m_inliers)[i]; }

  /*!
    \brief Returns the i^{th} item of the cluster.
  */
  const Item& operator[] (std::size_t i) const
  { return get (m_item_map, *(m_range->begin() + (*m_inliers)[i])); }

  /*!
    \brief Returns the bounding box of the cluster.
  */
  const CGAL::Bbox_3& bbox() const
  {
    if (m_bounding_box == CGAL::Bbox_3())
    {
      Point_idx_to_point_unary_function transform (m_range, m_item_map);
      m_bounding_box = CGAL::bbox_3 (boost::make_transform_iterator (m_inliers->begin(), transform),
                                     boost::make_transform_iterator (m_inliers->end(), transform));
    }

    return m_bounding_box;
  }

  /// @}

  /// \name Classification
  /// @{

  /*!
    \brief Returns the input classification value used for training.
  */
  int training() const { return m_training; }

  /*!
    \brief Returns a reference to the input classification value used for training.
  */
  int& training() { return m_training; }

  /*!
    \brief Returns the output classification value.
  */
  int label() const { return m_label; }

  /*!
    \brief Returns a reference to the output classification value.
  */
  int& label() { return m_label; }

  // @}
};

  /*!
    \ingroup PkgClassificationCluster

    \brief Given a set of cluster indices, segments the input `range`
    into `Cluster` objects.

    All items whose index value `idx` (accessed through `index_map`)
    is the same are stored in the same cluster at position `idx` in
    the `clusters` vector.

    \tparam ItemRange model of `ConstRange`. Its iterator type is
    `RandomAccessIterator`. Its value type depends on the data that is
    classified (for example, `CGAL::Point_3` or `CGAL::Triangle_3`).

    \tparam ItemMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `ItemRange` and value type
    is the type of item to classify (for example, `CGAL::Point_3`).

    \tparam IndexMap is a model of `ReadablePropertyMap` with value type `int`.

    \param range input range.
    \param item_map property map to access the input items.
    \param index_map property map that associates the index of an item
    in the input range to the index of a cluster (-1 if item is not
    assigned to a cluster).
    \param clusters container where generated `Cluster` objects are stored.
  */
template <typename ItemRange, typename ItemMap, typename IndexMap>
std::size_t create_clusters_from_indices (const ItemRange& range,
                                          ItemMap item_map,
                                          IndexMap index_map,
                                          std::vector<Cluster<ItemRange, ItemMap> >& clusters)
{
  std::size_t idx = 0;
  for (typename ItemRange::const_iterator it = range.begin(); it != range.end(); ++ it, ++ idx)
  {
    int c = int(get (index_map, idx));
    if (c == -1)
      continue;
    if (std::size_t(c) >= clusters.size())
    {
      clusters.reserve (c + 1);
      for (std::size_t i = clusters.size(); i <= std::size_t(c); ++ i)
        clusters.push_back (Cluster<ItemRange, ItemMap>(range, item_map));
    }
    clusters[std::size_t(c)].insert (idx);
  }

  return clusters.size();
}

} // namespace Classification

} // namespace CGAL


#endif // CGAL_CLASSIFICATION_CLUSTER_H
