// Copyright (c) 2011  INRIA Bordeaux Sud-Ouest (France), All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Gael Guennebaud

#ifndef CGAL_EIGEN_MATRIX_H
#define CGAL_EIGEN_MATRIX_H

#include <CGAL/basic.h> // include basic.h before testing #defines

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>

namespace CGAL {


/// The class Eigen_matrix
/// is a C++ wrapper around Eigen' matrix type SparseMatrix<>.
///
/// This kind of matrix can be either symmetric or not. Symmetric
/// matrices store only the lower triangle.
///
/// @heading Is Model for the Concepts: Model of the SparseLinearAlgebraTraits_d::Matrix concept.
///
/// @heading Parameters:
/// @param T Number type.

template<class T>
struct Eigen_matrix
{
// Public types
public:

   typedef Eigen::SparseMatrix<T> EigenType;
   typedef T NT;

// Public operations
public:

   /// Create a square matrix initialized with zeros.
   Eigen_matrix(int  dim,                  ///< Matrix dimension.
                bool is_symmetric = false) ///< Symmetric/hermitian?
     : m_dynamic(dim,dim), m_matrix(dim,dim)
   {
       CGAL_precondition(dim > 0);

       m_is_symmetric      = is_symmetric;
   }

   /// Create a rectangular matrix initialized with zeros.
   ///
   /// @commentheading Precondition: rows == columns if is_symmetric is true.
   Eigen_matrix(int  rows,                 ///< Number of rows.
                int  columns,              ///< Number of columns.
                bool is_symmetric = false) ///< Symmetric/hermitian?
     : m_dynamic(rows,columns), m_matrix(rows,columns)
   {
       CGAL_precondition(rows > 0);
       CGAL_precondition(columns > 0);
       if (m_is_symmetric) {
           CGAL_precondition(rows == columns);
       }

       m_is_symmetric      = is_symmetric;
   }

   /// Delete this object and the wrapped TAUCS matrix.
   ~Eigen_matrix()
   {
   }

   /// Return the matrix number of rows
   int row_dimension() const    { return m_matrix.rows(); }
   /// Return the matrix number of columns
   int column_dimension() const { return m_matrix.cols(); }


   /// Write access to a matrix coefficient: a_ij <- val.
   ///
   /// Optimizations:
   /// - For symmetric matrices, Eigen_matrix stores only the lower triangle
   ///   set_coef() does nothing if (i, j) belongs to the upper triangle.
   /// - Caller can optimize this call by setting 'new_coef' to true
   ///   if the coefficient does not already exist in the matrix.
   ///
   /// @commentheading Preconditions:
   /// - 0 <= i < row_dimension().
   /// - 0 <= j < column_dimension().
   void set_coef(int i, int j, T  val, bool new_coef = false)
   {
       CGAL_precondition(i < row_dimension());
       CGAL_precondition(j < column_dimension());

        if (m_is_symmetric && (j > i))
          return;

       m_is_uptodate = false;
       m_dynamic.coeffRef(i,j) = val;
   }



   const EigenType& eigen_object() const
   {
       if(!m_is_uptodate)
       {
         m_matrix = m_dynamic;
         m_is_uptodate = true;
       }
       return m_matrix;
   }

private:


   /// Eigen_matrix cannot be copied (yet)
   Eigen_matrix(const Eigen_matrix& rhs);
   Eigen_matrix& operator=(const Eigen_matrix& rhs);

// Fields
private:

   Eigen::DynamicSparseMatrix<T> m_dynamic;
   mutable EigenType m_matrix;
   mutable bool m_is_uptodate;

   // Symmetric/hermitian?
   bool m_is_symmetric;

}; // Taucs_matrix



/// The class Eigen_symmetric_matrix is a C++ wrapper
/// around a Eigen sparse matrix (type Eigen::SparseMatrix).
///
/// Symmetric matrices store only the lower triangle.
///
/// @heading Is Model for the Concepts: Model of the SparseLinearAlgebraTraits_d::Matrix concept.
///
/// @heading Parameters:
/// @param T Number type.

template<class T>
struct Eigen_symmetric_matrix
   : public Eigen_matrix<T>
{
// Public types
public:

   typedef T NT;

// Public operations
public:

   /// Create a square *symmetric* matrix initialized with zeros.
   Eigen_symmetric_matrix(int  dim)                  ///< Matrix dimension.
       : Eigen_matrix<T>(dim, true /* symmetric */)
   {
   }

   /// Create a square *symmetric* matrix initialized with zeros.
   ///
   /// @commentheading Precondition: rows == columns.
   Eigen_symmetric_matrix(int  rows,                 ///< Number of rows.
                          int  columns)              ///< Number of columns.
       : Eigen_matrix<T>(rows, columns, true /* symmetric */)
   {
   }
};

} //namespace CGAL

#endif // CGAL_EIGEN_MATRIX_H