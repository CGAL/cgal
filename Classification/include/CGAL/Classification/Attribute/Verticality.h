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

    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
    \tparam DiagonalizeTraits Solver used for matrix diagonalization.
  */
template <typename Kernel, typename RandomAccessIterator, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Verticality : public Attribute_base
{
  typedef Classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> verticality_attribute;
  
public:
  /*!
    \brief Constructs the attribute using local eigen analysis.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param eigen Class with precompute eigenvectors and eigenvalues
  */
  Verticality (RandomAccessIterator begin,
               RandomAccessIterator end,
               const Local_eigen_analysis& eigen)
  {
    this->weight() = 1.;
    typename Kernel::Vector_3 vertical (0., 0., 1.);

    for (std::size_t i = 0; i < (std::size_t)(end - begin); i++)
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

    \tparam VectorMap is a model of `ReadablePropertyMap` with value type `Vector_3<Kernel>`.
    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param normal_map Property map to access the normal vectors of the input points.
  */
  template <typename VectorMap>
  Verticality (const RandomAccessIterator& begin,
               const RandomAccessIterator& end,
               VectorMap normal_map)
  {
    this->weight() = 1.;
    typename Kernel::Vector_3 vertical (0., 0., 1.);

    for (std::size_t i = 0; i < (std::size_t)(end - begin); i++)
      {
        typename Kernel::Vector_3 normal = get(normal_map, begin[i]);
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

  virtual std::string id() { return "verticality"; }
  /// \endcond
};

} // namespace Attribute

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_ATTRIBUTE_VERTICALITY_H
