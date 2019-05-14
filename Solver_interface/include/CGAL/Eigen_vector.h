// Copyright (c) 2012  INRIA Bordeaux Sud-Ouest (France), All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Gael Guennebaud

#ifndef CGAL_EIGEN_VECTOR_H
#define CGAL_EIGEN_VECTOR_H

#include <CGAL/basic.h> // include basic.h before testing #defines

#include <Eigen/Core>

namespace CGAL {
/*!
\ingroup PkgSolver

The class `Eigen_vector` is a wrapper around \ref thirdpartyEigen "Eigen" vector
type <a href="http://eigen.tuxfamily.org/dox/classEigen_1_1Matrix.html"> </a>,
which is a simple array of numbers.

\cgalModels `SvdTraits::Vector`
\cgalModels `SparseLinearAlgebraTraits_d::Vector`.

\tparam T Number type.

\sa `CGAL::Eigen_solver_traits<T>`
\sa `CGAL::Eigen_sparse_matrix<T>`
\sa `CGAL::Eigen_sparse_symmetric_matrix<T>`
*/

template<class T>
class Eigen_vector
  : public Eigen::Matrix<T, Eigen::Dynamic, 1>
{
// Public types
public:
  /// \name Types
  /// @{
  typedef T                                      NT;

  /// The internal vector type from \ref thirdpartyEigen "Eigen".
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>    EigenType;
  /// @}

// Public operations
public:
  Eigen_vector<T>& operator=(const Eigen_vector<T>& other)
  {
    return static_cast<EigenType&>(*this) = other.eigen_object();
  }

  Eigen_vector<T>& operator=(const EigenType& other)
  {
    return static_cast<Eigen_vector<T>&>(static_cast<EigenType&>(*this) = other);
  }

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
