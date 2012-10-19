// Copyright (c) 2005  INRIA (France).
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

namespace SparseLinearAlgebraTraits_d {

/// \ingroup PkgSurfaceParameterizationConcepts
/// \cgalconcept
/// SparseLinearAlgebraTraits_d::Vector
/// is a concept of a vector that can be multiplied by a sparse matrix.
///
/// \refines LinearAlgebraTraits_d::Vector
/// \hasModel CGAL::Eigen_vector<T>
/// \hasModel OpenNL::FullVector<T> in OpenNL package
class Vector
{
// Public types
public:

    typedef Hidden_type NT;

// Public operations
public:

    /// Create a vector initialized with zeros.
    Vector (int rows);

    /// Copy constructor
    Vector(const Vector& toCopy);

    /// operator =()
    Vector& operator=(const Vector& toCopy);

    /// Return the vector's number of coefficients.
    int  dimension () const;

    /// Read/write access to a vector coefficient.
    ///
    /// \pre 0 <= row < dimension().
    NT  operator[] (int row) const;
    NT&  operator[] (int row);
};

/// \ingroup PkgSurfaceParameterizationConcepts
/// \cgalconcept
/// SparseLinearAlgebraTraits_d::Matrix
/// is a concept of a sparse matrix class.
///
/// \refines LinearAlgebraTraits_d::Matrix
/// \hasModel `CGAL::Eigen_sparse_matrix<T>`
/// \hasModel `CGAL::Eigen_sparse_symmetric_matrix<T>`
/// \hasModel `OpenNL::SparseMatrix<T>` in OpenNL package
///
/// \sa `SparseLinearAlgebraTraits_d`
/// \sa `SparseLinearAlgebraTraits_d::Vector`

class Matrix
{
// Public types
public:
    typedef Hidden_type NT;

// Public operations
public:
    /// Create a square matrix initialized with zeros.
    Matrix(int dimension);

    /// Create a rectangular matrix initialized with zeros.
    Matrix (int rows, int columns);

    /// Return the matrix number of rows.
    int  row_dimension () const;

    /// Return the matrix number of columns.
    int  column_dimension () const;

    /// Read access to a matrix coefficient.
    ///
    /// \pre 0 <= row < row_dimension().
    /// \pre 0 <= column < column_dimension().
    NT  get_coef (int row, int column) const;

    /// Write access to a matrix coefficient: a_ij <- a_ij + val.
    ///
    /// \pre 0 <= row < row_dimension().
    /// \pre 0 <= column < column_dimension().
    void add_coef(int row, int column, NT value);

    /// Write access to a matrix coefficient: a_ij <- val.
    ///
    /// Optimization:
    /// - Caller can optimize this call by setting `new_coef` to true
    ///   if the coefficient does not already exist in the matrix.
    ///
    /// \pre 0 <= i < row_dimension().
    /// \pre 0 <= j < column_dimension().
    void set_coef(int row, int column, NT value, bool new_coef = false);
};


}
