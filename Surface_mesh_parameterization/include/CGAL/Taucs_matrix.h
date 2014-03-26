// Copyright (c) 2005-2008  INRIA (France).
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
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy


#ifndef CGAL_TAUCS_MATRIX_H
#define CGAL_TAUCS_MATRIX_H

#include <CGAL/basic.h> // include basic.h before testing #defines

#ifdef CGAL_USE_TAUCS

#include <CGAL/Taucs_fix.h>
#include <CGAL/assertions.h>

#include <cstring>
#include <vector>

namespace CGAL {

/// @cond SKIP_IN_MANUAL
// Forward declaration
template<class T> struct Taucs_traits;


/// The class Taucs_matrix
/// is a C++ wrapper around TAUCS' matrix type taucs_ccs_matrix.
///
/// This kind of matrix can be either symmetric or not. Symmetric
/// matrices store only the lower triangle.
///
/// \cgalModels `SparseLinearAlgebraTraits_d::Matrix`
///
/// @param T Number type. Tested with T = taucs_single or taucs_double.
/// May also work with T = taucs_dcomplex and taucs_scomplex.

template<class T>
struct Taucs_matrix
{
// Public types
public:

    typedef T NT;

//__________________________________________________

// Private types
private:

   /**
    * A column of a Taucs_matrix. The column is
    * compressed, and stored in the form of
    * a vector of values + a vector of indices.
    */
    class Column
    {
    public:

        // Vector of values + vector of indices (linked)
        std::vector<T>   m_values;
        std::vector<int> m_indices;

    public:

        // Return the number of elements in the column
      int size() const    { return static_cast<int>(m_values.size()); }

        // return address of column{index}
        // (NULL if coefficient does not exist).
        const T* get_coef_addr(int index) const
        {
            // Search for element in vectors
            std::vector<int>::const_iterator        index_it;
            typename std::vector<T>::const_iterator value_it;
            for (index_it = m_indices.begin(), value_it = m_values.begin();
                 index_it != m_indices.end();
                 index_it++, value_it++)
            {
                if(*index_it == index)
                    return &*value_it;       // return address
            }

            // Element doesn't exist if we reach this point
            return NULL;
        }

        // column{index} <- column{index} + val
        void add_coef(int index, T val)
        {
            // Search for element in m_values[]
            T* coef_addr = (T*) get_coef_addr(index);
            if (coef_addr == NULL)
            {
              // if the coefficient doesn't exist yet
              m_indices.push_back(index);
              m_values.push_back(val);
            }
            else
            {
              // if the coefficient already exists
              *coef_addr += val;       // +=
            }
        }

        // column{index} <- val.
        void set_coef(int index, T val)
        {
            // Search for element in m_values[]
            T* coef_addr = (T*) get_coef_addr(index);
            if (coef_addr == NULL)
            {
              // if the coefficient doesn't exist yet
              m_indices.push_back(index);
              m_values.push_back(val);
            }
            else
            {
              // if the coefficient already exists
              *coef_addr = val;        // =
            }
        }

        // return column{index} (0 by default)
        T get_coef(int index) const
        {
            // Search for element in m_values[]
            const T* coef_addr = get_coef_addr(index);
            if (coef_addr == NULL)
              return 0; // if the coefficient doesn't exist
            else
              return *coef_addr;
        }
    }; // class Column

//__________________________________________________

// Public operations
public:

    /// Create a square matrix initialized with zeros.
    Taucs_matrix(int  dim,                  ///< Matrix dimension.
                 bool is_symmetric = false) ///< Symmetric/hermitian?
    {
        CGAL_precondition(dim > 0);

        m_row_dimension     = dim;
        m_column_dimension  = dim;
        m_is_symmetric      = is_symmetric;
        m_columns           = new Column[m_column_dimension];
        m_matrix            = NULL;
    }

    /// Create a rectangular matrix initialized with zeros.
    ///
    /// \pre rows == columns if is_symmetric is true.
    Taucs_matrix(int  rows,                 ///< Number of rows.
                 int  columns,              ///< Number of columns.
                 bool is_symmetric = false) ///< Symmetric/hermitian?
    {
        CGAL_precondition(rows > 0);
        CGAL_precondition(columns > 0);
        if (m_is_symmetric) {
            CGAL_precondition(rows == columns);
        }

        m_row_dimension     = rows;
        m_column_dimension  = columns;
        m_is_symmetric      = is_symmetric;
        m_columns           = new Column[m_column_dimension];
        m_matrix            = NULL;
    }

    /// Delete this object and the wrapped TAUCS matrix.
    ~Taucs_matrix()
    {
        // Delete the columns array
        delete[] m_columns;
        m_columns = NULL;

        // Delete the the wrapped TAUCS matrix
        if (m_matrix != NULL) {
            taucs_ccs_free(m_matrix);
            m_matrix = NULL;
        }
    }

    /// Return the matrix number of rows
    int row_dimension() const    { return m_row_dimension; }
    /// Return the matrix number of columns
    int column_dimension() const { return m_column_dimension; }

    /// Read access to a matrix coefficient.
    ///
    /// \pre 0 <= i < row_dimension().
    /// \pre 0 <= j < column_dimension().
    T  get_coef(int i, int j) const
    {
        CGAL_precondition(i < m_row_dimension);
        CGAL_precondition(j < m_column_dimension);

        // For symmetric matrices, we store only the lower triangle
        // => swap i and j if (i, j) belongs to the upper triangle
        if (m_is_symmetric && (j > i))
            std::swap(i, j);

        // Construct back the m_columns[] array after a call to get_taucs_matrix()
        if (m_columns == NULL)
          construct_back_columns();

        return m_columns[j].get_coef(i);
    }

    /// Write access to a matrix coefficient: a_ij <- val.
    ///
    /// Optimizations:
    /// - For symmetric matrices, Taucs_matrix stores only the lower triangle
    ///   set_coef() does nothing if (i, j) belongs to the upper triangle.
    /// - Caller can optimize this call by setting 'new_coef' to true
    ///   if the coefficient does not already exist in the matrix.
    ///
    /// \pre 0 <= i < row_dimension().
    /// \pre 0 <= j < column_dimension().
    void set_coef(int i, int j, T  val, bool new_coef = false)
    {
        CGAL_precondition(i < m_row_dimension);
        CGAL_precondition(j < m_column_dimension);

        if (m_is_symmetric && (j > i))
            return;

        // Construct back the m_columns[] array after a call to get_taucs_matrix()
        if (m_columns == NULL)
          construct_back_columns();

        // if caller knows that the coefficient doesn't exist yet
        if (new_coef)
        {
          m_columns[j].m_indices.push_back(i);
          m_columns[j].m_values.push_back(val);
        }
        else
        {
          m_columns[j].set_coef(i, val);
        }
    }

    /// Write access to a matrix coefficient: a_ij <- a_ij + val.
    ///
    /// Optimization:
    /// For symmetric matrices, Taucs_matrix stores only the lower triangle
    /// add_coef() does nothing if (i, j) belongs to the upper triangle.
    ///
    /// \pre 0 <= i < row_dimension().
    /// \pre 0 <= j < column_dimension().
    void add_coef(int i, int j, T val)
    {
        CGAL_precondition(i < m_row_dimension);
        CGAL_precondition(j < m_column_dimension);

        if (m_is_symmetric && (j > i))
            return;

        // Construct back the m_columns[] array after a call to get_taucs_matrix()
        if (m_columns == NULL)
          construct_back_columns();

        m_columns[j].add_coef(i, val);
    }

    /// Construct and return the TAUCS matrix wrapped by this object.
    /// The TAUCS matrix returned by this method is valid
    /// only until the next call to get_coef(), set_coef() or add_coef().
    //
    // Implementation note: this method deletes m_columns[] to save memory.
    const taucs_ccs_matrix* get_taucs_matrix() const
    {
        if (m_matrix == NULL)
        {
          CGAL_precondition(m_columns != NULL);

          // Convert matrix's T type to the corresponding TAUCS constant
          int flags = Taucs_traits<T>::TAUCS_FLAG;

          // We store only the lower triangle of symmetric matrices
          if (m_is_symmetric)
              flags |= TAUCS_TRIANGULAR | TAUCS_SYMMETRIC | TAUCS_LOWER;

          // Compute the number of non null elements in the matrix
          int nb_max_elements = 0;
          for (int col=0; col < m_column_dimension; col++)
              nb_max_elements += m_columns[col].size();

          // Allocate m_matrix
          m_matrix = taucs_ccs_create(m_row_dimension, m_column_dimension, nb_max_elements, flags);

          // Fill m_matrix. TAUCS matrix format is:
          // - rowind[] = array of non null elements of the matrix, ordered by columns
          // - values[] = array of row index of each element of rowind[]
          // - colptr[j] is the index of the first element of the column j (or where it
          //   should be if it doesn't exist) + the past-the-end index of the last column
          m_matrix->colptr[0] = 0;
          for (int col=0; col < m_column_dimension; col++)
          {
              int first_index = m_matrix->colptr[col]; // Index of 1st non null element of the column
              int nb_elements = m_columns[col].size(); // Number of non null elements of the column

              // Fast copy of column indices and values
	      std::memcpy(&m_matrix->rowind[first_index], &m_columns[col].m_indices[0], nb_elements*sizeof(int));
              T* taucs_values = (T*) m_matrix->values.v;
	      std::memcpy(&taucs_values[first_index], &m_columns[col].m_values[0],  nb_elements*sizeof(T));

              // Start of next column will be:
              m_matrix->colptr[col+1] = first_index + nb_elements;
          }

          // Delete m_columns[] to save memory.
          delete[] m_columns;
          m_columns = NULL;
        }

        CGAL_postcondition(m_matrix != NULL);
        CGAL_postcondition(m_columns == NULL);

        return m_matrix;
    }

private:

    // Construct back the m_columns[] array after a call to get_taucs_matrix(),
    // then delete m_matrix.
    void construct_back_columns() const
    {
        if (m_columns == NULL)
        {
          CGAL_precondition(m_matrix != NULL);

          // Allocate m_columns[].
          m_columns = new Column[m_column_dimension];

          // Fill m_columns[]. TAUCS matrix format is:
          // - rowind[] = array of non null elements of the matrix, ordered by columns
          // - values[] = array of row index of each element of rowind[]
          // - colptr[j] is the index of the first element of the column j (or where it
          //   should be if it doesn't exist) + the past-the-end index of the last column
          for (int col=0; col < m_column_dimension; col++)
          {
              int first_index = m_matrix->colptr[col]; // Index of 1st non null element of the column
              int nb_elements = m_matrix->colptr[col+1] - first_index;
                                                       // Number of non null elements of the column

              // Fast copy of column indices and values
              m_columns[col].m_indices.assign(&m_matrix->rowind[first_index],
                                              &m_matrix->rowind[first_index + nb_elements - 1]);
              T* taucs_values = (T*) m_matrix->values.v;
              m_columns[col].m_values.assign(&taucs_values[first_index],
                                             &taucs_values[first_index + nb_elements - 1]);
          }

          // Delete m_matrix
          taucs_ccs_free(m_matrix);
          m_matrix = NULL;
        }

        CGAL_postcondition(m_columns != NULL);
        CGAL_postcondition(m_matrix == NULL);
    }

    /// Taucs_matrix cannot be copied (yet)
    Taucs_matrix(const Taucs_matrix& rhs);
    Taucs_matrix& operator=(const Taucs_matrix& rhs);

// Fields
private:

    // Matrix dimensions
    int     m_row_dimension, m_column_dimension;

    // Symmetric/hermitian?
    bool    m_is_symmetric;

    // Matrix as a Columns array or a TAUCS matrix.
    // The matrix exists always as one of these kinds.
    mutable Column*           m_columns;
    mutable taucs_ccs_matrix* m_matrix;

}; // Taucs_matrix


/// The class Taucs_symmetric_matrix is a C++ wrapper
/// around a TAUCS *symmetric* matrix (type taucs_ccs_matrix).
///
/// Symmetric matrices store only the lower triangle.
///
/// \cgalModels `SparseLinearAlgebraTraits_d::Matrix`
///
/// @param T Number type. Tested with T = taucs_single or taucs_double.
/// May also work with T = taucs_dcomplex and taucs_scomplex.

template<class T>
struct Taucs_symmetric_matrix
    : public Taucs_matrix<T>
{
// Public types
public:

    typedef T NT;

// Public operations
public:

    /// Create a square *symmetric* matrix initialized with zeros.
    Taucs_symmetric_matrix(int  dim)                  ///< Matrix dimension.
        : Taucs_matrix<T>(dim, true /* symmetric */)
    {
    }

    /// Create a square *symmetric* matrix initialized with zeros.
    ///
    /// \pre rows == columns.
    Taucs_symmetric_matrix(int  rows,                 ///< Number of rows.
                           int  columns)              ///< Number of columns.
        : Taucs_matrix<T>(rows, columns, true /* symmetric */)
    {
    }
};


// Utility class for Taucs_matrix
// Convert matrix's T type to the corresponding TAUCS constant (called TAUCS_FLAG)
template<class T> struct Taucs_traits {};
template<> struct Taucs_traits<taucs_double> {
    enum { TAUCS_FLAG = TAUCS_DOUBLE };
};
template<> struct Taucs_traits<taucs_single>  {
    enum { TAUCS_FLAG = TAUCS_SINGLE };
};
template<> struct Taucs_traits<taucs_dcomplex> {
    enum { TAUCS_FLAG = TAUCS_DCOMPLEX };
};
template<> struct Taucs_traits<taucs_scomplex> {
    enum { TAUCS_FLAG = TAUCS_SCOMPLEX };
};
/// @endcond

} //namespace CGAL

#endif

#endif // CGAL_TAUCS_MATRIX_H
