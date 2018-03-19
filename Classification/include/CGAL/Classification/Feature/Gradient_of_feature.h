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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_GRADIENT_OF_FEATURE_H
#define CGAL_CLASSIFICATION_GRADIENT_OF_FEATURE_H

#include <CGAL/license/Classification.h>

#include <vector>

// Experimental feature, not used officially and not documented yet

/// \cond SKIP_IN_MANUAL

namespace CGAL {

namespace Classification {

namespace Feature {

  /*!
    \ingroup PkgClassificationFeatures

    TODO
  */
template <typename InputRange, typename ItemMap, typename NeighborQuery>
class Gradient_of_feature : public Feature_base
{
  const InputRange& m_input;
  ItemMap m_map;
  Feature_handle m_feature;
  boost::shared_ptr<NeighborQuery> m_query;
  
public:
  /*!
    TODO
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
  /// \cond SKIP_IN_MANUAL
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
  /// \endcond
};

} // namespace Feature

} // namespace Classification

} // namespace CGAL

/// \endcond

#endif // CGAL_CLASSIFICATION_GRADIENT_OF_FEATURE_H
