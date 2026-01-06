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

#ifndef CGAL_CLASSIFICATION_GRADIENT_OF_FEATURE_H
#define CGAL_CLASSIFICATION_GRADIENT_OF_FEATURE_H

#include <CGAL/license/Classification.h>

#include <vector>

// Experimental feature

namespace CGAL {

namespace Classification {

namespace Feature {

  /*!
    \ingroup PkgClassificationFeatures

    \brief Feature that computes the gradient of another feature.

    The gradient is estimated based on the local neighborhood of each item.
  */
template <typename InputRange, typename ItemMap, typename NeighborQuery>
class Gradient_of_feature : public Feature_base
{
  const InputRange& m_input;
  ItemMap m_map;
  Feature_handle m_feature;
  std::shared_ptr<NeighborQuery> m_query;

public:
  /*!
    \brief Constructs the feature.

    \param input the input range.
    \param map the property map.
    \param feature the feature whose gradient is computed.
    \param neighbor_query the neighbor query used to access the local neighborhood.
  */
  Gradient_of_feature (const InputRange& input,
                       ItemMap map,
                       Feature_handle feature,
                       const NeighborQuery& neighbor_query)
    : m_input (input), m_map (map), m_feature (feature), m_query (new NeighborQuery(neighbor_query))
  {
    std::ostringstream oss;
    oss << "gradient_of_" << feature->name();
    this->set_name (oss.str());
  }
  
  virtual float value (std::size_t pt_index)
  {
    std::vector<std::size_t> neighborhood;
    (*m_query) (get (m_map, *(m_input.begin()+pt_index)), std::back_inserter (neighborhood));

    if (neighborhood.empty())
      return 0.f;

    float mean = 0.f;

    for (std::size_t i = 0; i < neighborhood.size(); ++ i)
      mean += m_feature->value (neighborhood[i]);

    return (m_feature->value (pt_index) - mean / neighborhood.size());
  }
};

} // namespace Feature

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_GRADIENT_OF_FEATURE_H
