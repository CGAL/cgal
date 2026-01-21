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

#include <vector>
#include "accelerate/accelerate-swift.h"

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
  typedef SwiftAccelerate::Values Vec;
  /// @}

// Public operations
public:

  /// Constructs a null vector.
  Accelerate_vector()
    : m_vec(dimension), m_svec(Vec::init())
  {}

  /// Create a vector initialized with zeros.
  Accelerate_vector(std::size_t dimension)
    : m_vec(dimension), m_svec(Vec::init(dimension, dimension))
  {}

  void set(int i, NT v)
  {
    m_vec[i] = v;
  }

  /// Return the vector's number of coefficients.
  int dimension() const { return static_cast<int>(m_vec.size()); }

  NT& operator[](int row)
  {
    return m_vec[row];
  }

  void copy_back()
  {
    for(int i = 0; i < m_vec.size(); i++){
      m_vec[i] = m_svec.get(i);
    }
  }

  const NT& operator[](int row) const
  {
    return m_vec[row];
  }


  /*
  NT& operator[](int row)
  {
    return m_vec[row];
  }
  */
  /// Return a pointer to the data array of this vector.
  const Vec& data() const {
    for(int i = 0; i < m_vec.size(); i++){
      m_svec.set(i,m_vec[i]);
    }
    return m_svec;
  }

  Vec& data() {
    for(int i = 0; i < m_vec.size(); i++){
      m_svec.set(i,m_vec[i]);
    }
    return m_svec;
  }

private:
  std::vector<T> m_vec;
  mutable Vec  m_svec;
};

} // namespace CGAL

#endif // CGAL_ACCELERATE_VECTOR_H
