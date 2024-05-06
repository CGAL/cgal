
// Copyright (c) 2024  GeometryFactory SARL (France), All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_ACCELERATE_SOLVER_TRAIS_H
#define CGAL_ACCELERATE_SOLVER_TRAIS_H


#include <CGAL/Accelerate_vector.h>
#include <CGAL/Accelerate_sparse_matrix.h>

namespace CGAL {

/*!
  \ingroup PkgSolverInterfaceLS

  The class `Accelerate_solver_traits` provides an interface to the sparse solvers of
*/
template<class T>
class Accelerate_solver_traits
{
public:
  using NT = T;
  using Matrix = Accelerate_sparse_matrix<T>;
  using Vector = Accelerate_vector<T>;

  Accelerate_solver_traits()
    {}

  /// Solve the sparse linear system \f$ A \times X = B \f$.
  /// Return `true` on success. The solution is then \f$ (1/D) \times X \f$.
  ///
  /// \pre A.row_dimension() == B.dimension().
  /// \pre A.column_dimension() == X.dimension().
   bool linear_solver(const Matrix& A, const Vector& B, Vector& X, NT& D)
  {
    A.assemble_matrix();
    A.solve(B, X);
    X.copy_back();
    D = 1; // Accelerate does not support homogeneous coordinates
    return true;
  }

 };

} // namespace CGAL


#endif // ACCELERATE_SOLVER_TRAITS
