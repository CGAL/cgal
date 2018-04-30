// Copyright (c) 2018 GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_CLUSTER_H
#define CGAL_CLASSIFICATION_CLUSTER_H

#include <CGAL/license/Classification.h>


namespace CGAL {

namespace Classification {


template <typename ItemRange, typename ItemMap>
class Cluster
{
public:

  typedef typename ItemMap::value_type Item;
  
  struct Neighbor_query
  {
    template <typename OutputIterator>
    OutputIterator operator() (const Cluster& cluster, OutputIterator output) const
    {
      return std::copy (cluster.neighbors.begin(), cluster.neighbors.end(), output);
    }
  };
  
private:
  const ItemRange* m_range;
  ItemMap m_item_map;
  
  std::vector<std::size_t> m_inliers;
  mutable CGAL::Bbox_3 m_bounding_box;
  int m_training;
  int m_label;
  
public:
  std::vector<std::size_t> neighbors;
  
  Cluster (const ItemRange& range, ItemMap item_map)
    : m_range (&range), m_item_map (item_map)
    , m_training(-1), m_label(-1)
  { }

  void insert (std::size_t idx) { m_inliers.push_back (idx); }
  
  std::size_t size() const { return m_inliers.size(); }

  std::size_t index (std::size_t i) const { return m_inliers[i]; }
  const Item& operator[] (std::size_t i) const
  { return get (m_item_map, *(m_range->begin() + m_inliers[i])); }

  CGAL::Bbox_3 bbox() const
  {
    if (m_bounding_box == CGAL::Bbox_3())
    {
      m_bounding_box
        = CGAL::bbox_3 (boost::make_transform_iterator
                        (m_range->begin(),
                         CGAL::Property_map_to_unary_function<ItemMap>(m_item_map)),
                        boost::make_transform_iterator
                        (m_range->end(),
                         CGAL::Property_map_to_unary_function<ItemMap>(m_item_map)));

    }
    return m_bounding_box;
  }

  int training() const { return m_training; }
  int label() const { return m_label; }
  
  int& training() { return m_training; }
  int& label() { return m_label; }

};

template <typename ItemRange, typename ItemMap, typename IndexMap>
std::size_t create_clusters_from_indices (const ItemRange& range,
                                          ItemMap item_map,
                                          IndexMap index_map,
                                          std::vector<Cluster<ItemRange, ItemMap> >& clusters)
{
  std::size_t idx = 0;
  for (typename ItemRange::const_iterator it = range.begin(); it != range.end(); ++ it, ++ idx)
  {
    int c = int(get (index_map, *it));
    if (c == -1)
      continue;  
    if (std::size_t(c) >= clusters.size())
      clusters.resize (c + 1, Cluster<ItemRange, ItemMap>(range, item_map));
    clusters[std::size_t(c)].insert (idx);
  }
  
  return clusters.size();
}

} // namespace Classification

} // namespace CGAL


#endif // CGAL_CLASSIFICATION_CLUSTER_H
