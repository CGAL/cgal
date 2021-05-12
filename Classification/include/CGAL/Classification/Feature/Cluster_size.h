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

#ifndef CGAL_CLASSIFICATION_FEATURE_CLUSTER_SIZE_H
#define CGAL_CLASSIFICATION_FEATURE_CLUSTER_SIZE_H

#include <CGAL/license/Classification.h>

#include <vector>
#include <CGAL/Classification/Feature_base.h>

namespace CGAL {

namespace Classification {

namespace Feature {

  /*!
    \ingroup PkgClassificationCluster

    \brief %Feature that returns the size of each cluster.

    Its default name is "cluster_size".
  */
class Cluster_size : public CGAL::Classification::Feature_base
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
  Cluster_size (ClusterRange& clusters)
  {
    this->set_name ("cluster_size");
    m_values.reserve (clusters.size());

    for (std::size_t i = 0; i < clusters.size(); ++ i)
      m_values.push_back (float(clusters[i].size()));
  }

  /// \cond SKIP_IN_MANUAL
  virtual float value (std::size_t cluster_index)
  {
    return m_values[cluster_index];
  }
  /// \endcond
};

} // namespace Feature

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURE_CLUSTER_SIZE_H
