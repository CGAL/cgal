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


/// \ingroup PkgSurfaceParameterizationConcepts
/// \cgalconcept
/// The concept SparseLinearAlgebraTraits_d
/// is used to solve sparse linear systems "A*X = B".
/// \refines LinearAlgebraTraits_d
class SparseLinearAlgebraTraits_d
{
public:
    typedef xxx Matrix ;
    typedef xxx Vector ;
    typedef xxx NT;

// Public operations
public:
    /// Default constructor
    SparseLinearAlgebraTraits_d();

    /// Solve the sparse linear system "A*X = B".
    /// Return true on success. The solution is then (1/D) * X.
    ///
    /// \pre A.row_dimension() == B.dimension().
    /// \pre A.column_dimension() == X.dimension().
    bool linear_solver (const Matrix& A, const Vector& B, Vector& X, NT& D);
};

