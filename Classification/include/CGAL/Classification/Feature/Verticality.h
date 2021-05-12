// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_FEATURE_VERTICALITY_H
#define CGAL_CLASSIFICATION_FEATURE_VERTICALITY_H

#include <CGAL/license/Classification.h>

#include <vector>
#include <CGAL/Classification/Feature_base.h>
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

    Its default name is "verticality".

    \tparam GeomTraits model of \cgal Kernel.
  */
template <typename GeomTraits>
class Verticality : public Feature_base
{
  const typename GeomTraits::Vector_3 vertical;
  std::vector<compressed_float> verticality_feature;
  const Local_eigen_analysis* eigen;

public:
  /*!
    \brief constructs the feature using local eigen analysis.

    \tparam InputRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \param input input range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  template <typename InputRange>
  Verticality (const InputRange& input,
               const Local_eigen_analysis& eigen)
    : vertical (0., 0., 1.), eigen (&eigen)
  {
    CGAL_USE(input);
    this->set_name ("verticality");
  }


  /*!
    \brief constructs the feature using provided normals of points.

    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator` and its value type is the key type of
    `VectorMap`.
    \tparam VectorMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `GeomTraits::Vector_3`.
    \param input point range.
    \param normal_map property map to access the normal vectors of the input points.
  */
  template <typename PointRange, typename VectorMap>
  Verticality (const PointRange& input,
               VectorMap normal_map)
    : vertical (0., 0., 1.), eigen(nullptr)
  {
    this->set_name ("verticality");
    for (std::size_t i = 0; i < input.size(); i++)
    {
      typename GeomTraits::Vector_3 normal = get(normal_map, *(input.begin()+i));
      normal = normal / CGAL::sqrt (normal * normal);
      verticality_feature.push_back (compress_float(1.f - float(CGAL::abs(normal * vertical))));
    }
  }


  /// \cond SKIP_IN_MANUAL
  virtual float value (std::size_t pt_index)
  {
    if (eigen != nullptr)
    {
      typename GeomTraits::Vector_3 normal = eigen->normal_vector<GeomTraits>(pt_index);
      normal = normal / CGAL::sqrt (normal * normal);
      return (1.f - float(CGAL::abs(normal * vertical)));
    }
    else
      return decompress_float(verticality_feature[pt_index]);
  }
  /// \endcond
};

} // namespace Feature

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURE_VERTICALITY_H
