// Copyright (c) 2012  INRIA Bordeaux Sud-Ouest (France), All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

#ifndef CGAL_EIGEN_SOLVER_TRAITS_H
#define CGAL_EIGEN_SOLVER_TRAITS_H

#include <CGAL/basic.h> // include basic.h before testing #defines

#include <Eigen/Sparse>
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>

namespace CGAL {


namespace internal {
  template <class EigenSolver,class FT>
  struct Get_eigen_matrix{
    typedef Eigen_sparse_matrix<FT> type;
  };
  
  template <class FT,class EigenMatrix>
  struct Get_eigen_matrix< ::Eigen::ConjugateGradient<EigenMatrix>,FT>{
    typedef Eigen_sparse_symmetric_matrix<FT> type;
  };

  template <class FT,class EigenMatrix>
  struct Get_eigen_matrix< ::Eigen::SimplicialCholesky<EigenMatrix>,FT>{
    typedef Eigen_sparse_symmetric_matrix<FT> type;
  };
} //internal 
  
/// The class Eigen_solver_traits
/// is a generic traits class for solving asymmetric or symmetric positive definite (SPD)
/// sparse linear systems using one of the Eigen solvers.
/// The default solver is the iterative bi-congugate gradient stabilized solver
/// Eigen::BiCGSTAB for double.
///
/// @heading Is Model for the Concepts: Model of the SparseLinearAlgebraTraits_d concept.

template<class EigenSolverT = Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType> >
class Eigen_solver_traits
{
  typedef typename EigenSolverT::Scalar Scalar;
// Public types
public:
   typedef Scalar                                                       NT;
   typedef typename internal::Get_eigen_matrix<EigenSolverT,NT>::type   Matrix;
   typedef Eigen_vector<Scalar>                                         Vector;
   

// Public operations
public:

   Eigen_solver_traits()
   {
   }
   
   EigenSolverT& solver() { return m_solver; }

   /// Solve the sparse linear system "A*X = B".
   /// Return true on success. The solution is then (1/D) * X.
   ///
   /// @commentheading Preconditions:
   /// - A.row_dimension()    == B.dimension().
   /// - A.column_dimension() == X.dimension().
   bool linear_solver(const Matrix& A, const Vector& B, Vector& X, NT& D)
   {
      D = 1;          // Eigen does not support homogeneous coordinates

      m_solver.compute(A.eigen_object());
       
      if(m_solver.info() != Eigen::Success)
         return false;
         
      X = m_solver.solve(B);

      return m_solver.info() == Eigen::Success;
   }
protected:
  EigenSolverT m_solver;

};

} //namespace CGAL

#endif // CGAL_EIGEN_SOLVER_TRAITS_H
