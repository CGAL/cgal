// Copyright (c) 2016  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Simon Giraudot, Florent Lafarge

#ifndef CGAL_CLASSIFICATION_ATTRIBUTE_DISTANCE_TO_PLANE_H
#define CGAL_CLASSIFICATION_ATTRIBUTE_DISTANCE_TO_PLANE_H

#include <vector>

#include <CGAL/Classifier.h>

namespace CGAL {

namespace Classification {

namespace Attribute {

  /*!
    \ingroup PkgClassificationAttributes

    \brief Attribute based on local distance to a fitted plane.

    Characterizing a level of non-planarity can help identify noisy
    parts of the input such as vegetation. This attribute computes the
    distance of a point to a locally fitted plane.
    
    \tparam Kernel model of \cgal Kernel.
    \tparam Range range of items, model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `Range` and value type
    is `Point_3<Kernel>`.
    \tparam DiagonalizeTraits model of `DiagonalizeTraits` used
    for matrix diagonalization.
  */
template <typename Kernel, typename Range, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Distance_to_plane : public Attribute_base
{
  typedef Classification::Local_eigen_analysis<Kernel, Range,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;

  std::vector<double> distance_to_plane_attribute;
  
public:
  /*!
    \brief Constructs the attribute.

    \param input input range.
    \param point_map property map to access the input points.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  Distance_to_plane (const Range& input,
                     PointMap point_map,
                     const Local_eigen_analysis& eigen)
  {
    this->set_weight(1.);
    for(std::size_t i = 0; i < input.size(); i++)
      distance_to_plane_attribute.push_back
        (CGAL::sqrt (CGAL::squared_distance (get(point_map, input[i]), eigen.plane(i))));
    
    this->compute_mean_max (distance_to_plane_attribute, this->mean, this->max);
    //    max *= 2;
  }

  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return distance_to_plane_attribute[pt_index];
  }

  virtual std::string name() { return "distance_to_plane"; }
  /// \endcond
};

} // namespace Attribute

} // namespace Classification
  
} // namespace CGAL

#endif // CGAL_CLASSIFICATION_ATTRIBUTE_DISTANCE_TO_PLANE_H
