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
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Laurent Saboret, Pierre Alliez


#ifndef CGAL_TAUCS_MATRIX
#define CGAL_TAUCS_MATRIX

// Include TAUCS main header
#ifndef CGAL_INCLUDED_TAUCS_H
    #define CGAL_INCLUDED_TAUCS_H

    // GCC 3.3 complains if we include <complex.h> within "extern C {}"
    #if defined(__GNUC__)
        #include <complex.h>
    #endif

    // TAUCS is a C library
    extern "C" {
        #include <taucs.h>
    }

    // Avoid error with std::min() and std::max()
    #undef min
    #undef max

#endif

#include <cassert>

CGAL_BEGIN_NAMESPACE


// Forward declaration
template<class T> struct Taucs_number;


/// Struct Taucs_matrix
/// is a C++ wrapper around TAUCS' matrix type taucs_ccs_matrix.
///
/// This kind of matrix can be either symmetric or not. Symmetric
/// matrices store only the lower triangle.
///
/// TODO: reallocate the array of non null elements when it's full
///
/// Concept: Model of the SparseLinearAlgebraTraits_d::Matrix concept.

template<class T>       // Tested with T = taucs_single or taucs_double
                        // May also work with T = taucs_dcomplex and taucs_scomplex
struct Taucs_matrix
{
// Public types
public:

    typedef T NT;

// Public operations
public:

    /// Create a square matrix initialized with zeros
    Taucs_matrix(int  dim,                  ///< Matrix dimension
                 bool is_symmetric = false, ///< Symmetric/hermitian?
                 int  nb_max_elements = 0)  ///< Max number of non 0 elements in the
                                            ///< matrix (automatically computed if 0)
    {
        init(dim, dim, is_symmetric, nb_max_elements);
    }

    /// Create a rectangular matrix initialized with zeros
    Taucs_matrix(int  rows,                 ///< Matrix dimensions
                 int  columns,
                 bool is_symmetric = false, ///< Symmetric/hermitian?
                 int  nb_max_elements = 0)  ///< Max number of non 0 elements in the
                                            ///< matrix (automatically computed if 0)
    {
        init(rows, columns, is_symmetric, nb_max_elements);
    }

    ~Taucs_matrix() {
        /// Delete TAUCS matrix wrapped by this object
        taucs_ccs_free(m_matrix);
        m_matrix = NULL;
    }

    /// Return the matrix number of rows
    inline int row_dimension() const    {
        return m_matrix->m;
    }
    /// Return the matrix number of columns
    inline int column_dimension() const {
        return m_matrix->n;
    }

    /// Read access to 1 matrix coefficient
    ///
    /// Preconditions:
    /// - 0 <= i < row_dimension()
    /// - 0 <= j < column_dimension()
    T  get_coef (int i, int j) const
    {
        // For symmetric matrices, we store only the lower triangle
        // => swap i and j if (i,j) belongs to the upper triangle
        if (m_matrix->flags & TAUCS_SYMMETRIC)
        {
            assert(m_matrix->flags & TAUCS_LOWER);
            if (j > i)
                std::swap(i, j);
        }

        assert(i < row_dimension());
        assert(j < column_dimension());

        // Get pointer to matrix element (NULL if it doesn't exist)
        const T* element = find_element(i, j);

        return (element == NULL) ? 0 : (*element);
    }

    /// Write access to 1 matrix coefficient: aij <- val
    ///
    /// Optimization:
    /// For symmetric matrices, Taucs_matrix stores only the lower triangle
    /// set_coef() does nothing if (i,j) belongs to the upper triangle
    ///
    /// Preconditions:
    /// - 0 <= i < row_dimension()
    /// - 0 <= j < column_dimension()
    void set_coef(int i, int j, T  val)
    {
        if (m_matrix->flags & TAUCS_SYMMETRIC)
        {
            assert(m_matrix->flags & TAUCS_LOWER);
            if (j > i)
                return;
        }

        assert(i < row_dimension());
        assert(j < column_dimension());

        // Get pointer to matrix element. Create it if it doesn't exist.
        T* element = find_element(i, j, true);
        assert(element != NULL);

        *element = val;                     // =
    }

    /// Write access to 1 matrix coefficient: aij <- aij + val
    ///
    /// Optimization:
    /// For symmetric matrices, Taucs_matrix stores only the lower triangle
    /// add_coef() does nothing if (i,j) belongs to the upper triangle
    ///
    /// Preconditions:
    /// - 0 <= i < row_dimension()
    /// - 0 <= j < column_dimension()
    void add_coef(int i, int j, T val)
    {
        if (m_matrix->flags & TAUCS_SYMMETRIC)
        {
            assert(m_matrix->flags & TAUCS_LOWER);
            if (j > i)
                return;
        }

        assert(i < row_dimension());
        assert(j < column_dimension());

        // Get pointer to matrix element. Create it if it doesn't exist.
        T* element = find_element(i, j, true);
        assert(element != NULL);

        *element += val;                     // +=
    }

    /// Get TAUCS matrix wrapped by this object
    const taucs_ccs_matrix* get_taucs_matrix() const {
        return m_matrix;
    }
    taucs_ccs_matrix* get_taucs_matrix() {
        return m_matrix;
    }

// Private operations
private:

    /// Create a matrix initialized with zeros
    void init(int  rows,                    ///< Matrix dimensions
              int  columns,
              bool is_symmetric,            ///< Symmetric/hermitian?
              int  nb_max_elements)         ///< Max number of non 0 elements in the
                                            ///< matrix (automatically computed if 0)
    {
        assert(rows > 0);
        assert(columns > 0);

        // Convert matrix's T type to the corresponding TAUCS constant
        int flags = Taucs_number<T>::TAUCS_FLAG;

        // We store only the lower triangle of symmetric matrices
        if (is_symmetric)
        {
            assert(rows == columns);
            flags |= TAUCS_TRIANGULAR | TAUCS_SYMMETRIC | TAUCS_LOWER;
        }

        // Compute default max number of non null elements
        if (nb_max_elements <= 0)
        {
            // Pick a number larger than the average valence in a triangular mesh (6)
            int average_nb_elements_per_row_or_column = 16;
                                            // 16 is fine for parameterization package

            if (!is_symmetric)
                nb_max_elements = average_nb_elements_per_row_or_column * std::max(rows,columns);
            else
                nb_max_elements = (average_nb_elements_per_row_or_column/2+1) * std::max(rows,columns);
        }

        // Create TAUCS matrix wrapped by this object
        m_matrix = taucs_ccs_create(rows, columns, nb_max_elements, flags);

        // Init the current and max number of non null elements in m_matrix
        m_nb_elements = 0;
        m_nb_max_elements = nb_max_elements;

        // Init m_matrix's colptr[] array
        // Implementation note:
        // colptr[j] is the index of the first element of the column j (or where it
        // should be if it doesn't exist) + the past-the-end index of the last column
        for (int col=0; col <= columns; col++)
            m_matrix->colptr[col] = 0;
    }

    /// Read/write access to 1 matrix coefficient:
    /// Get a pointer to a matrix element. Optionaly create it.
    /// Return NULL if it doesn't exist (cannot happen if 'create' is true)
    ///
    /// Preconditions:
    /// - 0 <= i < row_dimension()
    /// - 0 <= j < column_dimension()
    /// - j <= i for symmetric matrices (we store only the lower triangle)
    T* find_element(int i, int j, bool create)
    {
        T* element = NULL;                          // returned value
        T* matrix_values = (T*)(m_matrix->values.v);// cast to actual array type

        assert(i < row_dimension());
        assert(j < column_dimension());
        if (m_matrix->flags & TAUCS_SYMMETRIC)
        {
            assert(m_matrix->flags & TAUCS_LOWER);
            assert(j <= i);
        }

        // Search for the i element of the j column's range in rowind[] and values[]
        // Implementation note:
        // colptr[j] is the index of the first element of the column j (or where it
        // should be if it doesn't exist) + the past-the-end index of the last column
        for (int idx = m_matrix->colptr[j]; idx <= m_matrix->colptr[j+1]; idx++)
        {
            // At the the (i,j) element's position (which may exist or not)
            if ((m_matrix->rowind[idx] >= i)    // on element i position of the column j
             || (idx == m_matrix->colptr[j+1])) // or on the 1st index of column j+1
            {
                // If the (i,j) element doesn't exist yet and 'create' is true,
                // shift the next elements and insert a new null element
                if (((m_matrix->rowind[idx] > i)    // too far in the column j
                 || (idx == m_matrix->colptr[j+1])) // or on the column j+1
                  && create )
                {
                    // TODO: reallocate m_matrix if it's full
                    assert(m_nb_elements < m_nb_max_elements);

                    // The number of elements to shift is:
                    int nb_elements_to_shift = m_matrix->colptr[column_dimension()] - idx;

                    // Shift values[] and insert a 0
                    if (nb_elements_to_shift > 0)
                        memmove(&matrix_values[idx+1], &matrix_values[idx], nb_elements_to_shift*sizeof(T));
                    matrix_values[idx] = 0;

                    // Shift rowind[] and add the i index
                    if (nb_elements_to_shift > 0)
                        memmove(&m_matrix->rowind[idx+1], &m_matrix->rowind[idx], nb_elements_to_shift*sizeof(int));
                    m_matrix->rowind[idx] = i;

                    // Shift colptr[]
                    for (int col=j+1; col <= column_dimension(); col++)
                        m_matrix->colptr[col]++;

                    m_nb_elements++;
                }

                // If the (i,j) element exists (now), return its address
                if ((m_matrix->rowind[idx] == i) && (idx < m_matrix->colptr[j+1]))
                    element = &matrix_values[idx];

                // exit for() statement
                break;
            }
        }

        return element;
    }

    /// Read access to 1 matrix coefficient:
    /// Get a pointer to a matrix element. Return NULL if it doesn't exist.
    ///
    /// Preconditions:
    /// - 0 <= i < row_dimension()
    /// - 0 <= j < column_dimension()
    /// - j <= i for symmetric matrices (we store only the lower triangle)
    const T* find_element(int i, int j) const
    {
        return ((Taucs_matrix<T>*)this)->find_element(i, j, false);
    }

    /// Taucs_matrix cannot be copied
    /// (for the moment, could be implemented if needed).
    Taucs_matrix(const Taucs_matrix& rhs);
    Taucs_matrix& operator=(const Taucs_matrix& rhs);

// Fields
private:

    /// The actual TAUCS matrix wrapped by this object
    taucs_ccs_matrix* m_matrix;

    /// Current and max number of non null elements in m_matrix
    int m_nb_elements;
    int m_nb_max_elements;
};


/// Struct Taucs_symmetric_matrix
/// is a C++ wrapper around a TAUCS *symmetric* matrix
/// (type taucs_ccs_matrix).
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

    /// Create a square SYMMETRIC matrix initialized with zeros
    Taucs_symmetric_matrix(int  dim,                  ///< Matrix dimension
                           int  nb_max_elements = 0)  ///< Max number of non 0 elements in the
                                                      ///< matrix (automatically computed if 0)
        : Taucs_matrix<T>(dim, true, nb_max_elements)
    {
    }

    /// Create a square SYMMETRIC matrix initialized with zeros
    Taucs_symmetric_matrix(int  rows,                 ///< Matrix dimensions
                           int  columns,
                           int  nb_max_elements = 0)  ///< Max number of non 0 elements in the
                                                      ///< matrix (automatically computed if 0)
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
