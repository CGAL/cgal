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

#include <CGAL/config.h> // include basic.h before testing #defines
#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4244)
#endif
#include <Eigen/Sparse>
#if EIGEN_VERSION_AT_LEAST(3, 1, 91)

#include <Eigen/SparseLU>
#endif


#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>
#include <boost/shared_ptr.hpp>

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
#if EIGEN_VERSION_AT_LEAST(3, 1, 91)
  template <class FT, class EigenMatrix, class EigenOrdering>
  struct Get_eigen_matrix< ::Eigen::SparseLU<EigenMatrix, EigenOrdering >, FT> {
    typedef Eigen_sparse_matrix<FT> type;
  };
#endif
} //internal 
  
/// The class Eigen_solver_traits
/// is a generic traits class for solving asymmetric or symmetric positive definite (SPD)
/// sparse linear systems using one of the Eigen solvers.
/// The default solver is the iterative bi-congugate gradient stabilized solver
/// <a href="http://eigen.tuxfamily.org/dox/classEigen_1_1BiCGSTAB.html">Eigen::BiCGSTAB</a> for double.
///
/// \cgalModels `SparseLinearAlgebraWithFactorTraits_d`.

template<class EigenSolverT = Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType> >
class Eigen_solver_traits
{
  typedef typename EigenSolverT::Scalar Scalar;
// Public types
public:
   typedef EigenSolverT Solver;
   typedef Scalar                                                       NT;
   typedef typename internal::Get_eigen_matrix<EigenSolverT,NT>::type   Matrix;
   typedef Eigen_vector<Scalar>                                         Vector;
   

// Public operations
public:

   Eigen_solver_traits():m_mat(NULL), m_solver_sptr(new EigenSolverT)
   {
   }
   
   EigenSolverT& solver() { return *m_solver_sptr; }

   /// Solve the sparse linear system "A*X = B".
   /// Return true on success. The solution is then (1/D) * X.
   ///
   /// @commentheading Preconditions:
   /// - A.row_dimension()    == B.dimension().
   /// - A.column_dimension() == X.dimension().
   bool linear_solver(const Matrix& A, const Vector& B, Vector& X, NT& D)
   {
      D = 1;          // Eigen does not support homogeneous coordinates

      m_solver_sptr->compute(A.eigen_object());
       
      if(m_solver_sptr->info() != Eigen::Success)
         return false;
         
      X = m_solver_sptr->solve(B);

      return m_solver_sptr->info() == Eigen::Success;
   }

  bool factor (const Matrix& A, NT& D)
  {
    D = 1;
    
    m_mat = &A.eigen_object();
    solver().compute(*m_mat);
    return solver().info() == Eigen::Success;
  }

  bool linear_solver(const Vector& B, Vector& X)
  {
    CGAL_precondition(m_mat!=NULL); //factor should have been called first
    X = solver().solve(B);
    return solver().info() == Eigen::Success;
  }

// Solving the normal equation "At*A*X = At*B".
// --
  bool normal_equation_factor(const Matrix& A)
  {
    typename Matrix::EigenType At = A.eigen_object().transpose();
    m_mat = &A.eigen_object();
    solver().compute(At * A.eigen_object());
    return solver().info() == Eigen::Success;
  }

  bool normal_equation_solver(const Vector& B, Vector& X)
  {
    CGAL_precondition(m_mat!=NULL); //non_symmetric_factor should have been called first
    typename Vector::EigenType AtB = m_mat->transpose() * B.eigen_object();
    X = solver().solve(AtB);
    return solver().info() == Eigen::Success;
  }

  bool normal_equation_solver(const Matrix& A, const Vector& B, Vector& X)
  {
    if (!normal_equation_factor(A)) return false;
    return normal_equation_solver(B, X);
  }
// --
protected:
  const typename Matrix::EigenType* m_mat;
  boost::shared_ptr<EigenSolverT> m_solver_sptr;

};

//specialization of the solver for BiCGSTAB as for surface parameterization, the 
//intializer should be a vector of one's (this was the case in 3.1-alpha but not in the official 3.1).
template<>
class Eigen_solver_traits< Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType> >
{
  typedef Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType> EigenSolverT;
  typedef EigenSolverT::Scalar Scalar;
// Public types
public:
   typedef EigenSolverT Solver;
   typedef Scalar                                                       NT;
   typedef internal::Get_eigen_matrix<EigenSolverT,NT>::type   Matrix;
   typedef Eigen_vector<Scalar>                                         Vector;
   

// Public operations
public:

   Eigen_solver_traits(): m_solver_sptr(new EigenSolverT)
   {
   }
   
   EigenSolverT& solver() { return *m_solver_sptr; }

   /// Solve the sparse linear system "A*X = B".
   /// Return true on success. The solution is then (1/D) * X.
   ///
   /// @commentheading Preconditions:
   /// - A.row_dimension()    == B.dimension().
   /// - A.column_dimension() == X.dimension().
   bool linear_solver(const Matrix& A, const Vector& B, Vector& X, NT& D)
   {
      D = 1;          // Eigen does not support homogeneous coordinates

      m_solver_sptr->compute(A.eigen_object());
       
      if(m_solver_sptr->info() != Eigen::Success)
         return false;
      
      X.setOnes(B.rows());
      X = m_solver_sptr->solveWithGuess(B,X);

      return m_solver_sptr->info() == Eigen::Success;
   }
protected:
  boost::shared_ptr<EigenSolverT> m_solver_sptr;

};

} //namespace CGAL

#endif // CGAL_EIGEN_SOLVER_TRAITS_H
