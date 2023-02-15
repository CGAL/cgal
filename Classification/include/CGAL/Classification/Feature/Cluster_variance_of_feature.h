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

#ifndef CGAL_CLASSIFICATION_FEATURE_CLUSTER_VARIANCE_FEATURE_H
#define CGAL_CLASSIFICATION_FEATURE_CLUSTER_VARIANCE_FEATURE_H

#include <CGAL/license/Classification.h>

#include <vector>
#include <sstream>

#include <CGAL/Classification/Feature_base.h>

namespace CGAL {

namespace Classification {

namespace Feature {

  /*!
    \ingroup PkgClassificationCluster

    \brief %Feature that computes the variance values of an itemwise
    feature over the respective items of clusters.

    Its default name is "variance_" + the name of the itemwise feature.
  */
class Cluster_variance_of_feature : public CGAL::Classification::Feature_base
{
  std::vector<float> m_values;

public:

  /*!
    \brief constructs the feature.

    \tparam ClusterRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator` and its value type is the key type of
    `Cluster`.

    \param clusters input range.
    \param itemwise_feature feature that takes values on the range of
    items from which `clusters` is a subset.
    \param mean_feature `Cluster_mean_of_feature` computed on
    `itemwise_feature`.
  */
  template <typename ClusterRange>
  Cluster_variance_of_feature (ClusterRange& clusters,
                               Feature_handle itemwise_feature,
                               Feature_handle mean_feature)
  {
    std::ostringstream oss;
    oss << "variance_" << itemwise_feature->name();
    this->set_name (oss.str());

    m_values.reserve (clusters.size());
    for (std::size_t i = 0; i < clusters.size(); ++ i)
    {
      double mean = double (mean_feature->value(i));
      double variance = 0.;

      for (std::size_t j = 0; j < clusters[i].size(); ++ j)
      {
        double v = double (itemwise_feature->value (clusters[i].index(j)));
        variance += (v - mean) * (v - mean);
      }
      variance /= clusters[i].size();
      m_values.push_back (float(variance));
    }
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

#endif // CGAL_CLASSIFICATION_FEATURE_CLUSTER_VARIANCE_FEATURE_H
