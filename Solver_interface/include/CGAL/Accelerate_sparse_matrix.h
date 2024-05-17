// Copyright (c) 2024  GeometryFactory SARL (France), All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_SOLVER_INTERFACE_ACCELERATE_SPARSE_MATRIX_H
#define CGAL_SOLVER_INTERFACE_ACCELERATE_SPARSE_MATRIX_H

#include <CGAL/basic.h> // include basic.h before testing #defines
#include <CGAL/Accelerate_vector.h>
#include "accelerate/accelerate-swift.h"

#include <map>

namespace CGAL {

/*!
\ingroup PkgSolverInterfaceLS

The class `Accelerate_sparse_matrix` is a wrapper around the `SparseMatrix_Double` matrix type
<a href=https://developer.apple.com/documentation/accelerate/sparsematrix_double">`Accelerate::SparseMatrix`_Double</a>
that represents general matrices, be they symmetric or not.

\cgalModels{SparseLinearAlgebraTraits_d::Matrix}

\tparam T Number type.

\sa `CGAL::Accelerate_vector<T>`
\sa `CGAL::Accelerate_matrix<T>`
\sa `CGAL::Accelerate_sparse_symmetric_matrix<T>`
*/
template<class T>
class Accelerate_sparse_matrix
{
  // Public types
public:
  /// \name Types
  /// @{

  /// The internal matrix type from \ref thirdpartyAccelerate "Accelerate".
  typedef SwiftAccelerate::Matrix Matrix;

  typedef T                           NT;
  /// @}

#if 0
  Accelerate_sparse_matrix(const AccelerateType& et)
    : m_is_already_built(true), m_has_been_changed(false), m_matrix(et), m_is_symmetric(false)
  {}
#endif

  // Public operations
public:
  Accelerate_sparse_matrix()
    : m_matrix(Matrix::init()), m_is_already_built(false), m_has_been_changed(false)
  { }

  /// Create a square matrix initialized with zeros.
  Accelerate_sparse_matrix(int dim,            ///< Matrix dimension.
                           bool is_symmetric = false)  ///< Symmetric/hermitian?
    :  m_matrix(Matrix::init()), m_rows(dim), m_columns(dim), m_is_already_built(false), m_has_been_changed(false),  m_is_symmetric(is_symmetric)
  {
    CGAL_precondition(dim > 0);
  }


  /// Create a rectangular matrix initialized with zeros.
  ///
  /// \pre rows == columns if `is_symmetric` is true.
  Accelerate_sparse_matrix(int rows,          ///< Number of rows.
                           int columns,       ///< Number of columns.
                           bool is_symmetric = false) ///< Symmetric/hermitian?
    :  m_matrix(Matrix::init()), m_rows(rows), m_columns(columns), m_is_already_built(false), m_has_been_changed(false), m_is_symmetric(is_symmetric)
  {
    CGAL_precondition(rows > 0);
    CGAL_precondition(columns > 0);
    if(m_is_symmetric)
    {
      CGAL_precondition(rows == columns);
    }

    // reserve memory for a regular 3D grid
  }

  void swap(Accelerate_sparse_matrix& other)
  {
    std::swap(m_rows, other.m_rows);
    std::swap(m_columns, other.m_columns);
    std::swap(m_is_already_built, other.m_is_already_built);
    std::swap(m_has_been_changed, other.m_has_been_changed);
    std::swap(m_is_symmetric, other.m_is_symmetric);
    //m_matrix.swap(other.m_matrix);
    m_triplets.swap(other.m_triplets);
  }

  /// Delete this object and the wrapped matrix.
  ~Accelerate_sparse_matrix() = default;

  /// Return the matrix number of rows
  int row_dimension() const { return m_rows; }

  /// Return the matrix number of columns
  int column_dimension() const { return m_columns; }

  /// Write access to a matrix coefficient: a_ij <- val.
  ///
  /// Users can optimize calls to this function by setting 'new_coef' to `true`
  /// if the coefficient does not already exist in the matrix.
  ///
  /// \warning For symmetric matrices, `Accelerate_sparse_matrix` only stores the lower triangle
  ///   and `set_coef()` does nothing if (i, j) belongs to the upper triangle.
  ///
  /// \pre 0 <= i < row_dimension().
  /// \pre 0 <= j < column_dimension().
  void set_coef(std::size_t i_, std::size_t j_, T val, bool new_coef = false)
  {
    int i = static_cast<int>(i_);
    int j = static_cast<int>(j_);
    CGAL_precondition(i < row_dimension());
    CGAL_precondition(j < column_dimension());

    if (m_is_symmetric && (j > i))
      return;

    m_triplets[std::make_pair(i,j)] = val;
    if(m_is_already_built){
      m_has_been_changed = true;
    }
  }

  /// Write access to a matrix coefficient: a_ij <- a_ij + val.
  ///
  /// \warning For symmetric matrices, Accelerate_sparse_matrix only stores the lower triangle
  /// `add_coef()` does nothing if (i, j) belongs to the upper triangle.
  ///
  /// \pre 0 <= i < row_dimension().
  /// \pre 0 <= j < column_dimension().
  void add_coef(std::size_t i_, std::size_t j_, T val)
  {
    int i = static_cast<int>(i_);
    int j = static_cast<int>(j_);
    if(m_is_symmetric && (j > i))
      return;

    m_triplets[std::make_pair(i,j)] += val;
    if(m_is_already_built){
      m_has_been_changed = true;
    }
  }

  /// Read access to a matrix coefficient.
  ///
  /// \pre 0 <= i < row_dimension().
  /// \pre 0 <= j < column_dimension().
  NT get_coef (std::size_t i_, std::size_t j_) const
  {
    int i = static_cast<int>(i_);
    int j = static_cast<int>(j_);
    CGAL_precondition(i < row_dimension());
    CGAL_precondition(j < column_dimension());

    if(m_is_symmetric && j > i){
      std::swap(i, j);
    }
    auto it = m_triplets.find(std::make_pair(i,j));
    if(it != m_triplets.end()){
      return it->second;
    }
    return 0;
  }

  /// \cond SKIP_IN_MANUAL
  void assemble_matrix() const
  {
    SwiftAccelerate::Indices rows = SwiftAccelerate::Indices::init();
    SwiftAccelerate::Indices columns = SwiftAccelerate::Indices::init();
    SwiftAccelerate::Values values = SwiftAccelerate::Values::init();

    for(const auto& triplet : m_triplets){
      int i = triplet.first.first;
      int j = triplet.first.second;
      const NT& val = triplet.second;
      rows.push_back(i);
      columns.push_back(j);
      values.push_back(val);
      if(m_is_symmetric && (i != j)){
        rows.push_back(j);
        columns.push_back(i);
        values.push_back(val);
      }
    }
    m_matrix.initialize(m_rows, rows, columns, values);
    m_is_already_built = true;
    m_has_been_changed = false;
  }
  /// \endcond

public:

  void solve(const Accelerate_vector<T>& B, Accelerate_vector<T>& X) const
  {
    Accelerate_vector<T>& ncB = const_cast<Accelerate_vector<T>&>(B);
    Accelerate_vector<T>& ncX = const_cast<Accelerate_vector<T>&>(X);
    m_matrix.solve(ncB.data(), ncX.data());
  }
  
  /// \cond SKIP_IN_MANUAL
  friend Accelerate_sparse_matrix
  operator*(const T& c, const Accelerate_sparse_matrix& M)
  {
    std::cout << "todo: implement operator*(scalar, matrix)" << std::endl;
    assert(false);
    return M;
  }

  friend Accelerate_sparse_matrix
  operator+(const Accelerate_sparse_matrix& M0, const Accelerate_sparse_matrix& M1)
  {
    std::cout << "todo: implement operator+( matrix, matrix)" << std::endl;
    assert(false);
    return M0;
  }
  /// \endcond

  // Fields
private:
  mutable bool m_is_already_built;
  mutable bool m_has_been_changed;
  bool m_is_symmetric;
  int m_rows, m_columns;
  mutable std::map<std::pair<int,int>,NT> m_triplets;

  mutable Matrix m_matrix;

}; // Accelerate_sparse_matrix




} //namespace CGAL

#endif // CGAL_SOLVER_INTERFACE_ACCELERATE_SPARSE_MATRIX_H
