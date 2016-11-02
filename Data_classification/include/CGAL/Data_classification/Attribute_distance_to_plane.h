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

#ifndef CGAL_DATA_CLASSIFICATION_ATTRIBUTE_DISTANCE_TO_PLANE_H
#define CGAL_DATA_CLASSIFICATION_ATTRIBUTE_DISTANCE_TO_PLANE_H

#include <vector>

#include <CGAL/Point_set_classification.h>

namespace CGAL {

namespace Data_classification {

  /*!
    \ingroup PkgDataClassification

    \brief Attribute based on local distance to a fitted plane.

    Characterizing a level of non-planarity can help identify noisy
    parts of the input such as vegetation. This attribute computes the
    distance of a point to a locally fitted plane.
    
    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
    \tparam DiagonalizeTraits Solver used for matrix diagonalization.
  */
template <typename Kernel, typename RandomAccessIterator, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Attribute_distance_to_plane : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointMap, DiagonalizeTraits> Local_eigen_analysis;

  std::vector<double> distance_to_plane_attribute;
  
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_map Property map to access the input points
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Attribute_distance_to_plane (RandomAccessIterator begin,
                               RandomAccessIterator end,
                               PointMap point_map,
                               const Local_eigen_analysis& eigen)
  {
    this->weight = 1.;
    for(std::size_t i = 0; i < (std::size_t)(end - begin); i++)
      distance_to_plane_attribute.push_back
        (CGAL::sqrt (CGAL::squared_distance (get(point_map, begin[i]), eigen.plane(i))));
    
    this->compute_mean_max (distance_to_plane_attribute, this->mean, this->max);
    //    max *= 2;
  }

  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return distance_to_plane_attribute[pt_index];
  }

  virtual std::string id() { return "distance_to_plane"; }
  /// \endcond
};


} // namespace Data_classification
  
} // namespace CGAL

#endif // CGAL_DATA_CLASSIFICATION_ATTRIBUTE_DISTANCE_TO_PLANE_H
