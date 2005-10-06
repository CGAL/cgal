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


/// Concept SparseLinearAlgebraTraits_d::Matrix
/// is a concept of a sparse matrix class.
///
/// Sub-concept: This is a sub-concept of LinearAlgebraTraits_d::Matrix.
///
/// Models:
/// - Taucs_matrix
/// - OpenNL::SparseMatrix

class Matrix
{
// Public types
public:
    typedef xxx NT;

// Public operations
public:
    /// Create a square matrix initialized with zeros
    Matrix(int dimension);

    /// Create a rectangular matrix initialized with zeros
    Matrix (int rows, int columns);

    /// Return the matrix number of rows
    int  row_dimension () const;

    /// Return the matrix number of columns
    int  column_dimension () const;

    /// Read access to 1 matrix coefficient.
    ///
    /// Preconditions:
    /// - 0 <= row < row_dimension().
    /// - 0 <= column < column_dimension().
    NT  get_coef (int row, int column) const;

    /// Write access to 1 matrix coefficient: aij <- aij + val.
    ///
    /// Preconditions:
    /// - 0 <= row < row_dimension().
    /// - 0 <= column < column_dimension().
    void add_coef(int row, int column, NT value);

    /// Write access to 1 matrix coefficient.
    ///
    /// Preconditions:
    /// - 0 <= row < row_dimension().
    /// - 0 <= column < column_dimension().
    void  set_coef (int row, int column, NT value);
};

