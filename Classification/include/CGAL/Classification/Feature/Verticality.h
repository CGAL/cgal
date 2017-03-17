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
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_FEATURE_VERTICALITY_H
#define CGAL_CLASSIFICATION_FEATURE_VERTICALITY_H

#include <vector>

#include <CGAL/Classification/Local_eigen_analysis.h>

namespace CGAL {

namespace Classification {

namespace Feature {

  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on local verticality. The orientation of the
    local tangent plane of the considered point can be useful to
    discriminate facades from the ground. This feature can use
    normal vectors if available or eigen analysis that estimates
    tangent planes from local neighborhoods.

    \tparam Geom_traits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `Geom_traits::Point_3`.
    \tparam DiagonalizeTraits model of `DiagonalizeTraits` used
    for matrix diagonalization.
  */
template <typename Geom_traits, typename PointRange, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<float,3> >
class Verticality : public Feature_base
{
  typedef Classification::Local_eigen_analysis<Geom_traits, PointRange,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;
  
  const typename Geom_traits::Vector_3 vertical;
  std::vector<float> verticality_feature;
  const Local_eigen_analysis* eigen;
  
public:
  /*!
    \brief Constructs the feature using local eigen analysis.

    \param input input range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
  Verticality (const PointRange& input,
               const Local_eigen_analysis& eigen)
    : vertical (0., 0., 1.), eigen(NULL)
  {
    this->set_name ("verticality");

    for (std::size_t i = 0; i < input.size(); i++)
      {
        typename Geom_traits::Vector_3 normal = eigen.normal_vector(i);
        normal = normal / CGAL::sqrt (normal * normal);
        verticality_feature.push_back (1. - CGAL::abs(normal * vertical));
      }
  }
#else
  Verticality (const PointRange&,
               const Local_eigen_analysis& eigen)
    : vertical (0., 0., 1.), eigen (&eigen)
  {
    this->set_name ("verticality");
  }
#endif

  /*!
    \brief Constructs the feature using provided normals of points.

    \tparam VectorMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `Geom_traits::Vector_3`.
    \param input input range.
    \param normal_map property map to access the normal vectors of the input points.
  */
  template <typename VectorMap>
  Verticality (const PointRange& input,
               VectorMap normal_map)
    : vertical (0., 0., 1.), eigen(NULL)
  {
    this->set_name ("verticality");
    for (std::size_t i = 0; i < input.size(); i++)
      {
        typename Geom_traits::Vector_3 normal = get(normal_map, *(input.begin()+i));
        normal = normal / CGAL::sqrt (normal * normal);
        verticality_feature.push_back (1. - std::fabs(normal * vertical));
      }
  }


  /// \cond SKIP_IN_MANUAL
  virtual float value (std::size_t pt_index)
  {
#ifndef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
    if (eigen != NULL)
      {
        typename Geom_traits::Vector_3 normal = eigen->normal_vector(pt_index);
        normal = normal / CGAL::sqrt (normal * normal);
        return (1. - CGAL::abs(normal * vertical));
      }
    else
#endif
      return verticality_feature[pt_index];
  }
  /// \endcond
};

} // namespace Feature

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURE_VERTICALITY_H
