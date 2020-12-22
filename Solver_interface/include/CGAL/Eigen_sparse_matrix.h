// Copyright (c) 2012  INRIA Bordeaux Sud-Ouest (France), All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Gael Guennebaud

#ifndef CGAL_SOLVER_INTERFACE_EIGEN_SPARSE_MATRIX_H
#define CGAL_SOLVER_INTERFACE_EIGEN_SPARSE_MATRIX_H

#include <CGAL/basic.h> // include basic.h before testing #defines

#include <Eigen/Sparse>

namespace CGAL {

/*!
\ingroup PkgSolverInterfaceRef

The class `Eigen_sparse_matrix` is a wrapper around `Eigen` matrix type
<a href="http://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html">`Eigen::SparseMatrix`</a>
that represents general matrices, be they symmetric or not.

\cgalModels `SparseLinearAlgebraTraits_d::Matrix`

\tparam T Number type.

\sa `CGAL::Eigen_vector<T>`
\sa `CGAL::Eigen_matrix<T>`
\sa `CGAL::Eigen_sparse_symmetric_matrix<T>`
*/
template<class T>
struct Eigen_sparse_matrix
{
  // Public types
public:
  /// \name Types
  /// @{

  /// The internal matrix type from \ref thirdpartyEigen "Eigen".
  typedef Eigen::SparseMatrix<T>      EigenType;

  typedef T                           NT;
  /// @}

  Eigen_sparse_matrix(const EigenType& et)
    : m_is_already_built(true), m_matrix(et), m_is_symmetric(false)
  {}

  // Public operations
public:
  Eigen_sparse_matrix() : m_is_already_built(false) { }

  /// Create a square matrix initialized with zeros.
  Eigen_sparse_matrix(std::size_t dim,            ///< Matrix dimension.
                      bool is_symmetric = false)  ///< Symmetric/hermitian?
    :
      m_is_already_built(false),
      m_matrix(static_cast<int>(dim), static_cast<int>(dim))
  {
    CGAL_precondition(dim > 0);

    m_is_symmetric = is_symmetric;
    m_triplets.reserve(dim); // reserve memory for a regular 3D grid
  }

  /// Create a square matrix initialized with zeros.
  Eigen_sparse_matrix(int dim,                    ///< Matrix dimension.
                      bool is_symmetric = false)  ///< Symmetric/hermitian?
    : m_is_already_built(false),
      m_matrix(dim, dim)
  {
    CGAL_precondition(dim > 0);

    m_is_symmetric = is_symmetric;
    // reserve memory for a regular 3D grid
    m_triplets.reserve(dim);
  }

  /// Create a rectangular matrix initialized with zeros.
  ///
  /// \pre rows == columns if `is_symmetric` is true.
  Eigen_sparse_matrix(std::size_t rows,          ///< Number of rows.
                      std::size_t columns,       ///< Number of columns.
                      bool is_symmetric = false) ///< Symmetric/hermitian?
    : m_is_already_built(false),
      m_matrix(static_cast<int>(rows), static_cast<int>(columns))
  {
    CGAL_precondition(rows > 0);
    CGAL_precondition(columns > 0);
    if(m_is_symmetric)
    {
      CGAL_precondition(rows == columns);
    }

    m_is_symmetric = is_symmetric;
    // reserve memory for a regular 3D grid
    m_triplets.reserve(rows);
  }

  void swap(Eigen_sparse_matrix& other)
  {
    std::swap(m_is_already_built, other.m_is_already_built);
    std::swap(m_is_symmetric, other.m_is_symmetric);
    m_matrix.swap(other.m_matrix);
    m_triplets.swap(other.m_triplets);
  }

  /// Delete this object and the wrapped matrix.
  ~Eigen_sparse_matrix() { }

  /// Create a rectangular matrix initialized with zeros.
  ///
  /// \pre rows == columns if `is_symmetric` is true.
  Eigen_sparse_matrix(int rows,                  ///< Number of rows.
                      int columns,               ///< Number of columns.
                      bool is_symmetric = false) ///< Symmetric/hermitian?
    : m_is_already_built(false),
      m_matrix(rows,columns)
  {
    CGAL_precondition(rows > 0);
    CGAL_precondition(columns > 0);
    if(is_symmetric)
    {
      CGAL_precondition(rows == columns);
    }

    m_is_symmetric = is_symmetric;
    // reserve memory for a regular 3D grid
    m_triplets.reserve(rows);
  }

  /// Return the matrix number of rows
  int row_dimension() const { return static_cast<int>(m_matrix.rows()); }
  /// Return the matrix number of columns
  int column_dimension() const { return static_cast<int>(m_matrix.cols()); }

  /// Write access to a matrix coefficient: a_ij <- val.
  ///
  /// Users can optimize calls to this function by setting 'new_coef' to `true`
  /// if the coefficient does not already exist in the matrix.
  ///
  /// \warning For symmetric matrices, `Eigen_sparse_matrix` only stores the lower triangle
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

    if(m_is_already_built)
    {
      m_matrix.coeffRef(i,j) = val;
    }
    else
    {
      if(new_coef == false)
      {
        assemble_matrix();
        m_matrix.coeffRef(i,j) = val;
      }
      else
      {
        m_triplets.push_back(Triplet(i,j,val));
      }
    }
  }

  /// Write access to a matrix coefficient: a_ij <- a_ij + val.
  ///
  /// \warning For symmetric matrices, Eigen_sparse_matrix only stores the lower triangle
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

    if(m_is_already_built){
      CGAL_precondition(i < row_dimension());
      CGAL_precondition(j < column_dimension());
      m_matrix.coeffRef(i,j) += val;
    }else{
      m_triplets.push_back(Triplet(i,j,val));
    }
  }

  /// Read access to a matrix coefficient.
  ///
  /// \warning Complexity:
  /// - O(log(n)) if the matrix is already built.
  /// - O(n) if the matrix is not built.
  /// `n` being the number of entries in the matrix.
  ///
  /// \pre 0 <= i < row_dimension().
  /// \pre 0 <= j < column_dimension().
  NT get_coef (std::size_t i_, std::size_t j_) const
  {
    int i = static_cast<int>(i_);
    int j = static_cast<int>(j_);
    CGAL_precondition(i < row_dimension());
    CGAL_precondition(j < column_dimension());

    if(m_is_symmetric && j > i)
      std::swap(i, j);

    if (m_is_already_built)
      return m_matrix.coeffRef(i,j);
    else
    {
      NT val = 0;
      for(std::size_t t=0; t<m_triplets.size(); ++t)
      {
        if(m_triplets[t].col() == j &&
           m_triplets[t].row() == i)
          val += m_triplets[t].value();
      }
      return val;
    }
  }

  /// \cond SKIP_IN_MANUAL
  void assemble_matrix() const
  {
    m_matrix.setFromTriplets(m_triplets.begin(), m_triplets.end());
    m_is_already_built = true;
    m_triplets.clear(); // the matrix is built and will not be rebuilt
  }
  /// \endcond

  /// Return the internal matrix, with type `EigenType`.
  const EigenType& eigen_object() const
  {
    if(!m_is_already_built)
      assemble_matrix();

    // turns the matrix into compressed mode:
    //  -> release some memory
    //  -> required for some external solvers
    m_matrix.makeCompressed();
    return m_matrix;
  }

  /// Return the internal matrix, with type `EigenType`.
  EigenType& eigen_object()
  {
    if(!m_is_already_built)
      assemble_matrix();

    // turns the matrix into compressed mode:
    //  -> release some memory
    //  -> required for some external solvers
    m_matrix.makeCompressed();
    return m_matrix;
  }

public:
  /// \cond SKIP_IN_MANUAL
  friend Eigen_sparse_matrix
  operator*(const T& c, const Eigen_sparse_matrix& M)
  {
    return Eigen_sparse_matrix(c* M.eigen_object());
  }

  friend Eigen_sparse_matrix
  operator+(const Eigen_sparse_matrix& M0, const Eigen_sparse_matrix& M1)
  {
    return Eigen_sparse_matrix(M0.eigen_object()+ M1.eigen_object());
  }
  /// \endcond

  // Fields
private:
  mutable bool m_is_already_built;

  typedef Eigen::Triplet<T,int> Triplet;
  mutable std::vector<Triplet> m_triplets;

  mutable EigenType m_matrix;

  // Symmetric/hermitian?
  bool m_is_symmetric;
}; // Eigen_sparse_matrix


/*!
\ingroup PkgSolverInterfaceRef

The class `Eigen_sparse_symmetric_matrix` is a wrapper around `Eigen` matrix type
<a href="http://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html">`Eigen::SparseMatrix` </a>

Since the matrix is symmetric, only the lower triangle part is stored.

\cgalModels `SparseLinearAlgebraTraits_d::Matrix`

\tparam T Number type.

\sa `CGAL::Eigen_vector<T>`
\sa `CGAL::Eigen_sparse_matrix<T>`
*/
template<class T>
struct Eigen_sparse_symmetric_matrix
  : public Eigen_sparse_matrix<T>
{
  /// Create a square *symmetric* matrix initialized with zeros.
  Eigen_sparse_symmetric_matrix(int dim)                  ///< Matrix dimension.
    : Eigen_sparse_matrix<T>(dim, true /* symmetric */)
  { }

  /// Create a square *symmetric* matrix initialized with zeros.
  ///
  /// \pre rows == columns.
  Eigen_sparse_symmetric_matrix(int rows,                 ///< Number of rows.
                                int columns)              ///< Number of columns.
    : Eigen_sparse_matrix<T>(rows, columns, true /* symmetric */)
  { }
};

} //namespace CGAL

#endif // CGAL_SOLVER_INTERFACE_EIGEN_SPARSE_MATRIX_H
