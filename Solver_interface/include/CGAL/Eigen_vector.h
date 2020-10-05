// Copyright (c) 2012  INRIA Bordeaux Sud-Ouest (France), All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Gael Guennebaud

#ifndef CGAL_EIGEN_VECTOR_H
#define CGAL_EIGEN_VECTOR_H

#include <CGAL/basic.h> // include basic.h before testing #defines

#include <Eigen/Core>

namespace CGAL {
/*!
\ingroup PkgSolverInterfaceRef

The class `Eigen_vector` is a wrapper around `Eigen`
<a href="http://eigen.tuxfamily.org/dox/classEigen_1_1Matrix.html">vector type</a>,
which is a simple array of numbers.

\cgalModels `SvdTraits::Vector`
\cgalModels `SparseLinearAlgebraTraits_d::Vector`.

\tparam T Number type.

\sa `CGAL::Eigen_solver_traits<T>`
\sa `CGAL::Eigen_sparse_matrix<T>`
\sa `CGAL::Eigen_sparse_symmetric_matrix<T>`
*/

template<class T, int D = ::Eigen::Dynamic>
class Eigen_vector
  : public ::Eigen::Matrix<T, D, 1>
{
// Public types
public:
  /// \name Types
  /// @{
  typedef T                                      NT;

  /// The internal vector type from \ref thirdpartyEigen "Eigen".
  typedef ::Eigen::Matrix<T, D, 1>               EigenType;
  /// @}

// Public operations
public:
  Eigen_vector& operator=(const Eigen_vector& other)
  {
    return static_cast<EigenType&>(*this) = other.eigen_object();
  }

  Eigen_vector& operator=(const EigenType& other)
  {
    return static_cast<Eigen_vector&>(static_cast<EigenType&>(*this) = other);
  }

  /// Constructs a null vector.
  Eigen_vector() : EigenType() {}

  /// Create a vector initialized with zeros.
  Eigen_vector(std::size_t dimension)
    : EigenType(static_cast<int>(dimension))
  {
    this->setZero();
  }

  /// Copy constructor.
  Eigen_vector(const Eigen_vector& toCopy) : EigenType(toCopy) { }

  ~Eigen_vector() { }

  /// Return the vector's number of coefficients.
  int dimension() const { return static_cast<int>(this->size()); }

  /// Return the internal vector wrapped by this object.
  const EigenType& eigen_object() const { return *this; }

  /// Return the internal vector wrapped by this object.
  EigenType& eigen_object() { return *this; }

  /// Write access to a vector coefficient: `a_i` <- `value`.
  void set(std::size_t i, NT value)
  {
    this->operator[](static_cast<int>(i)) = value;
  }

  /// Return a pointer to the data array of this vector.
  NT* vector() { return this->data(); }
};

} // namespace CGAL

#endif // CGAL_EIGEN_VECTOR_H
