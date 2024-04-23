// Copyright (c) 2024  GeometryFactory SARL (France), All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_ACCELERATE_VECTOR_H
#define CGAL_ACCELERATE_VECTOR_H


namespace CGAL {
/*!
\ingroup PkgSolverInterfaceLS

The class `Accelerate_vector` is a vector of numbers.

\cgalModels{SvdTraits::Vector,SparseLinearAlgebraTraits_d::Vector}

\tparam T Number type. Must be `double` or `float`

\sa `CGAL::Accelerate_solver_traits<T>`
\sa `CGAL::Accelerate_sparse_matrix<T>`
*/

template<class T>
class Accelerate_vector

{
// Public types
public:
  /// \name Types
  /// @{
  typedef T                                      NT;

  /// @}

// Public operations
public:

  /// Constructs a null vector.
  Accelerate_vector() = default;

  /// Create a vector initialized with zeros.
  Accelerate_vector(std::size_t dimension)
    : m_vec(dimension, NT(0))
  {}


  /// Return the vector's number of coefficients.
  int dimension() const { return static_cast<int>(m_vec.size()); }

  NT operator[](int row) const
  {
    return m_vec[row];
  }


  NT& operator[](int row)
  {
    return m_vec[row];
  }

  /// Return a pointer to the data array of this vector.
  const NT* data() { return &m_vec[0]; }

private:
  std::vector<NT> m_vec;
};

} // namespace CGAL

#endif // CGAL_ACCELERATE_VECTOR_H
