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

#ifndef CGAL_CLASSIFICATION_FEATURE_CLUSTER_VERTICAL_EXTENT_H
#define CGAL_CLASSIFICATION_FEATURE_CLUSTER_VERTICAL_EXTENT_H

#include <CGAL/license/Classification.h>

#include <vector>

#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Bbox_3.h>

namespace CGAL {

namespace Classification {

namespace Feature {

  /*!
    \ingroup PkgClassificationCluster

    \brief %Feature that returns the length of the smallest interval
    on the `Z` axis that contains all the items of a cluster.

    Its default name is "cluster_vertical_extent".
  */
class Cluster_vertical_extent : public CGAL::Classification::Feature_base
{
  std::vector<float> m_values;

public:

  /*!
    \brief constructs the feature.

    \tparam ClusterRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator` and its value type is the key type of
    `Cluster`.

    \param clusters input range.
  */
  template <typename ClusterRange>
  Cluster_vertical_extent (const ClusterRange& clusters)
  {
    typedef typename ClusterRange::const_iterator::value_type::Item Item;

    this->set_name ("cluster_vertical_extent");

    m_values.reserve (clusters.size());
    for (std::size_t i = 0; i < clusters.size(); ++ i)
    {
      float min_z = (std::numeric_limits<float>::max)();
      float max_z = -(std::numeric_limits<float>::min)();

      for (std::size_t j = 0; j < clusters[i].size(); ++ j)
      {
        const Item& item = clusters[i][j];
        const CGAL::Bbox_3& bbox = item.bbox();
        min_z = (std::min) (float(bbox.zmin()), min_z);
        max_z = (std::max) (float(bbox.zmax()), max_z);
      }
      m_values.push_back ((max_z - min_z));
    }
  }

  /// \cond SKIP_IN_MANUAL
  virtual float value (std::size_t cluster_index) { return m_values[cluster_index]; }
  /// \endcond

};

} // namespace Feature

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURE_CLUSTER_VERTICAL_EXTENT_H
