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


/// SparseLinearAlgebraTraits_d::Matrix
/// is a concept of a sparse matrix class.
///
/// @heading Has Models:
/// - Taucs_matrix<T>
/// - Taucs_symmetric_matrix<T>
/// - OpenNL::SparseMatrix<T> in OpenNL package

class Matrix : public LinearAlgebraTraits_d::Matrix
{
// Public types
public:
    typedef xxx NT;

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
    /// @commentheading Preconditions:
    /// - 0 <= row < row_dimension().
    /// - 0 <= column < column_dimension().
    NT  get_coef (int row, int column) const;

    /// Write access to a matrix coefficient: a_ij <- a_ij + val.
    ///
    /// @commentheading Preconditions:
    /// - 0 <= row < row_dimension().
    /// - 0 <= column < column_dimension().
    void add_coef(int row, int column, NT value);

    /// Write access to a matrix coefficient: a_ij <- val.
    ///
    /// Optimization:
    /// - Caller can optimize this call by setting 'new_coef' to true
    ///   if the coefficient does not already exist in the matrix.
    ///
    /// @commentheading Preconditions:
    /// - 0 <= i < row_dimension().
    /// - 0 <= j < column_dimension().
    void set_coef(int row, int column, NT value, bool new_coef = false);
};


} // namespace SparseLinearAlgebraTraits_d
