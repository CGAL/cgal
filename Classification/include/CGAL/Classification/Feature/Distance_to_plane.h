// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
// Copyright (c) 2017 GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot, Florent Lafarge

#ifndef CGAL_CLASSIFICATION_FEATURE_DISTANCE_TO_PLANE_H
#define CGAL_CLASSIFICATION_FEATURE_DISTANCE_TO_PLANE_H

#include <CGAL/license/Classification.h>
#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Classification/Local_eigen_analysis.h>
#include <vector>

namespace CGAL {

namespace Classification {

namespace Feature {

  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on local distance to a fitted
    plane. Characterizing a level of non-planarity can help to
    identify noisy parts of the input such as vegetation. This
    feature computes the distance of a point to a locally fitted
    plane.

    Its default name is "distance_to_plane".
    
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator` and its value type is the key type of
    `PointMap`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `CGAL::Point_3`.
  */
template <typename PointRange, typename PointMap>
class Distance_to_plane : public Feature_base
{

  typedef typename CGAL::Kernel_traits<typename PointMap::value_type>::Kernel Kernel;
  
#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
  std::vector<float> distance_to_plane_feature;
#else
  const PointRange& input;
  PointMap point_map;
  const Local_eigen_analysis& eigen;
#endif
  
public:
  /*!
    \brief Constructs the feature.

    \param input point range.
    \param point_map property map to access the input points.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  Distance_to_plane (const PointRange& input,
                     PointMap point_map,
                     const Local_eigen_analysis& eigen)
#ifndef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
    : input(input), point_map(point_map), eigen(eigen)
#endif
  {
    this->set_name ("distance_to_plane");
#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES    
    for(std::size_t i = 0; i < input.size(); i++)
      distance_to_plane_feature.push_back
        (CGAL::sqrt (CGAL::squared_distance (get(point_map, *(input.begin()+i)), eigen.plane<Kernel>(i))));
#endif
  }

  /// \cond SKIP_IN_MANUAL
  virtual float value (std::size_t pt_index)
  {
#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
    return distance_to_plane_feature[pt_index];
#else
    return float(CGAL::sqrt (CGAL::squared_distance
                             (get(point_map, *(input.begin()+pt_index)),
                              eigen.plane<Kernel>(pt_index))));
#endif
  }
  /// \endcond
};

} // namespace Feature

} // namespace Classification
  
} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURE_DISTANCE_TO_PLANE_H
