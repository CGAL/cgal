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
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_ATTRIBUTE_VERTICALITY_H
#define CGAL_CLASSIFICATION_ATTRIBUTE_VERTICALITY_H

#include <vector>

#include <CGAL/Classification/Local_eigen_analysis.h>

namespace CGAL {

namespace Classification {

namespace Attribute {

  /*!
    \ingroup PkgClassificationAttributes

    \brief Attribute based on local verticality.

    The orientation of the best fitting plane of a local neighborhood
    of the considered point can be useful to discriminate facades from
    the ground.

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
class Verticality : public Attribute_base
{
  typedef Classification::Local_eigen_analysis<Kernel, Range,
                                               PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> verticality_attribute;
  
public:
  /*!
    \brief Constructs the attribute using local eigen analysis.

    \param input input range.
    \param eigen class with precomputed eigenvectors and eigenvalues.
  */
  Verticality (const Range& input,
               const Local_eigen_analysis& eigen)
  {
    this->set_weight(1.);
    typename Kernel::Vector_3 vertical (0., 0., 1.);

    for (std::size_t i = 0; i < input.size(); i++)
      {
        typename Kernel::Vector_3 normal = eigen.normal_vector(i);
        normal = normal / CGAL::sqrt (normal * normal);
        verticality_attribute.push_back (1. - CGAL::abs(normal * vertical));
      }
    
    this->compute_mean_max (verticality_attribute, this->mean, this->max);
    //    max *= 2;
  }

  /*!
    \brief Constructs the attribute using provided normals of points.

    \tparam VectorMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `Range` and value type
    is `Vector_3<Kernel>`.
    \param input input range.
    \param normal_map property map to access the normal vectors of the input points.
  */
  template <typename VectorMap>
  Verticality (const Range& input,
               VectorMap normal_map)
  {
    this->set_weight(1.);
    typename Kernel::Vector_3 vertical (0., 0., 1.);

    for (std::size_t i = 0; i < input.size(); i++)
      {
        typename Kernel::Vector_3 normal = get(normal_map, *(input.begin()+i));
        normal = normal / CGAL::sqrt (normal * normal);
        verticality_attribute.push_back (1. - std::fabs(normal * vertical));
      }
    
    this->compute_mean_max (verticality_attribute, this->mean, this->max);
    //    max *= 2;
  }


  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return verticality_attribute[pt_index];
  }

  virtual std::string name() { return "verticality"; }
  /// \endcond
};

} // namespace Attribute

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_ATTRIBUTE_VERTICALITY_H
