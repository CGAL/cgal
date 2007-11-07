// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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


#ifndef CGAL_TAUCS_MATRIX
#define CGAL_TAUCS_MATRIX

#include <CGAL/basic.h>
#include <CGAL/Taucs_fix.h>


CGAL_BEGIN_NAMESPACE


// Forward declaration
template<class T> struct Taucs_number;


/// The class Taucs_matrix
/// is a C++ wrapper around TAUCS' matrix type taucs_ccs_matrix.
///
/// This kind of matrix can be either symmetric or not. Symmetric
/// matrices store only the lower triangle.
///
/// Concept: Model of the SparseLinearAlgebraTraits_d::Matrix concept.
///
template<class T>       // Tested with T = taucs_single or taucs_double
                        // May also work with T = taucs_dcomplex and taucs_scomplex
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
        int dimension() const    { return m_values.size(); }

        // column{index} <- column{index} + val
        void add_coef(int index, T val)
        {
            // Search for element in vectors
            std::vector<int>::iterator          index_it;
            typename std::vector<T>::iterator   value_it;
            for (index_it = m_indices.begin(), value_it = m_values.begin();
                 index_it != m_indices.end();
                 index_it++, value_it++)
            {
                if(*index_it == index) {
                    *value_it += val;       // +=
                    return;
                }
            }

            // Element doesn't exist yet if we reach this point
            m_indices.push_back(index);
            m_values.push_back(val);
        }

        // column{index} <- val
        void set_coef(int index, T val)
        {
            // Search for element in vectors
            std::vector<int>::iterator          index_it;
            typename std::vector<T>::iterator   value_it;
            for (index_it = m_indices.begin(), value_it = m_values.begin();
                 index_it != m_indices.end();
                 index_it++, value_it++)
            {
                if(*index_it == index) {
                    *value_it = val;        // =
                    return;
                }
            }

            // Element doesn't exist yet if we reach this point
            m_indices.push_back(index);
            m_values.push_back(val);
        }

        // return column{index} (0 by default)
        T get_coef(int index) const
        {
            // Search for element in vectors
            std::vector<int>::const_iterator        index_it;
            typename std::vector<T>::const_iterator value_it;
            for (index_it = m_indices.begin(), value_it = m_values.begin();
                 index_it != m_indices.end();
                 index_it++, value_it++)
            {
                if(*index_it == index)
                    return *value_it;       // return value
            }

            // Element doesn't exist yet if we reach this point
            return 0;
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
        m_columns           = new Column[m_column_dimension];
        m_is_symmetric      = is_symmetric;
        m_matrix            = NULL;
    }

    /// Create a rectangular matrix initialized with zeros.
    Taucs_matrix(int  rows,                 ///< Matrix dimensions.
                 int  columns,
                 bool is_symmetric = false) ///< Symmetric/hermitian?
    {
        CGAL_precondition(rows > 0);
        CGAL_precondition(columns > 0);
        if (m_is_symmetric) {
            CGAL_precondition(rows == columns);
        }

        m_row_dimension     = rows;
        m_column_dimension  = columns;
        m_columns           = new Column[m_column_dimension];
        m_is_symmetric      = is_symmetric;
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
    /// Preconditions:
    /// - 0 <= i < row_dimension().
    /// - 0 <= j < column_dimension().
    T  get_coef(int i, int j) const
    {
        // For symmetric matrices, we store only the lower triangle
        // => swap i and j if (i, j) belongs to the upper triangle
        if (m_is_symmetric && (j > i))
            std::swap(i, j);

        CGAL_precondition(i < m_row_dimension);
        CGAL_precondition(j < m_column_dimension);
        return m_columns[j].get_coef(i);
    }

    /// Write access to a matrix coefficient: a_ij <- val.
    ///
    /// Optimization:
    /// For symmetric matrices, Taucs_matrix stores only the lower triangle
    /// set_coef() does nothing if (i, j) belongs to the upper triangle.
    ///
    /// Preconditions:
    /// - 0 <= i < row_dimension().
    /// - 0 <= j < column_dimension().
    void set_coef(int i, int j, T  val)
    {
        if (m_is_symmetric && (j > i))
            return;

        CGAL_precondition(i < m_row_dimension);
        CGAL_precondition(j < m_column_dimension);
        m_columns[j].set_coef(i, val);
    }

    /// Write access to a matrix coefficient: a_ij <- a_ij + val.
    ///
    /// Optimization:
    /// For symmetric matrices, Taucs_matrix stores only the lower triangle
    /// add_coef() does nothing if (i, j) belongs to the upper triangle.
    ///
    /// Preconditions:
    /// - 0 <= i < row_dimension().
    /// - 0 <= j < column_dimension().
    void add_coef(int i, int j, T val)
    {
        if (m_is_symmetric && (j > i))
            return;

        CGAL_precondition(i < m_row_dimension);
        CGAL_precondition(j < m_column_dimension);
        m_columns[j].add_coef(i, val);
    }

    /// Construct and return the TAUCS matrix wrapped by this object.
    /// Note: the TAUCS matrix returned by this method is valid
    ///       only until the next call to set_coef(), add_coef() or get_taucs_matrix().
    const taucs_ccs_matrix* get_taucs_matrix() const
    {
        if (m_matrix != NULL) {
            taucs_ccs_free(m_matrix);
            m_matrix = NULL;
        }

        // Convert matrix's T type to the corresponding TAUCS constant
        int flags = Taucs_number<T>::TAUCS_FLAG;

        // We store only the lower triangle of symmetric matrices
        if (m_is_symmetric)
            flags |= TAUCS_TRIANGULAR | TAUCS_SYMMETRIC | TAUCS_LOWER;

        // Compute the number of non null elements in the matrix
        int nb_max_elements = 0;
        for (int col=0; col < m_column_dimension; col++)
            nb_max_elements += m_columns[col].dimension();

        // Create the TAUCS matrix wrapped by this object
        m_matrix = taucs_ccs_create(m_row_dimension, m_column_dimension, nb_max_elements, flags);

        // Fill m_matrix's colptr[], rowind[] and values[] arrays
        // Implementation note:
        // - rowind[] = array of non null elements of the matrix, ordered by columns
        // - values[] = array of row index of each element of rowind[]
        // - colptr[j] is the index of the first element of the column j (or where it
        //   should be if it doesn't exist) + the past-the-end index of the last column
        m_matrix->colptr[0] = 0;
        for (int col=0; col < m_column_dimension; col++)
        {
            // Number of non null elements of the column
            int nb_elements = m_columns[col].dimension();

            // Fast copy of column indices and values
            memcpy(&m_matrix->rowind[m_matrix->colptr[col]], &m_columns[col].m_indices[0], nb_elements*sizeof(int));
            T* taucs_values = (T*) m_matrix->values.v;
            memcpy(&taucs_values[m_matrix->colptr[col]], &m_columns[col].m_values[0],  nb_elements*sizeof(T));

            // Start of next column will be:
            m_matrix->colptr[col+1] = m_matrix->colptr[col] + nb_elements;
        }

        return m_matrix;
    }

private:

    /// Taucs_matrix cannot be copied (yet)
    Taucs_matrix(const Taucs_matrix& rhs);
    Taucs_matrix& operator=(const Taucs_matrix& rhs);

// Fields
private:

    // Matrix dimensions
    int     m_row_dimension, m_column_dimension;

    // Columns array
    Column* m_columns;

    // Symmetric/hermitian?
    bool    m_is_symmetric;

    /// The actual TAUCS matrix wrapped by this object.
    // This is in fact a COPY of the columns array
    mutable taucs_ccs_matrix* m_matrix;

}; // Taucs_matrix


/// The class Taucs_symmetric_matrix is a C++ wrapper
/// around a TAUCS *symmetric* matrix (type taucs_ccs_matrix).
///
/// Symmetric matrices store only the lower triangle.
///
/// Concept: Model of the SparseLinearAlgebraTraits_d::Matrix concept.

template<class T>       // Tested with T = taucs_single or taucs_double
                        // May also work with T = taucs_dcomplex and taucs_scomplex
struct Taucs_symmetric_matrix
    : public Taucs_matrix<T>
{
// Public types
public:

    typedef T NT;

// Public operations
public:

    /// Create a square SYMMETRIC matrix initialized with zeros.
    /// The max number of non 0 elements in the matrix is automatically computed.
    Taucs_symmetric_matrix(int  dim)                  ///< Matrix dimension.
        : Taucs_matrix<T>(dim, true)
    {
    }

    /// Create a square SYMMETRIC matrix initialized with zeros.
    Taucs_symmetric_matrix(int  rows,                 ///< Matrix dimensions.
                           int  columns,
                           int  nb_max_elements = 0)  ///< Max number of non 0 elements in the
                                                      ///< matrix (automatically computed if 0).
        : Taucs_matrix<T>(rows, columns, true, nb_max_elements)
    {
    }
};


// Utility class to Taucs_matrix
// Convert matrix's T type to the corresponding TAUCS constant (called TAUCS_FLAG)
template<class T> struct Taucs_number {};
template<> struct Taucs_number<taucs_double> {
    enum { TAUCS_FLAG = TAUCS_DOUBLE };
};
template<> struct Taucs_number<taucs_single>  {
    enum { TAUCS_FLAG = TAUCS_SINGLE };
};
template<> struct Taucs_number<taucs_dcomplex> {
    enum { TAUCS_FLAG = TAUCS_DCOMPLEX };
};
template<> struct Taucs_number<taucs_scomplex> {
    enum { TAUCS_FLAG = TAUCS_SCOMPLEX };
};


CGAL_END_NAMESPACE

#endif // CGAL_TAUCS_MATRIX
