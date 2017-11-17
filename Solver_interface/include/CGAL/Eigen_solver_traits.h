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
// SPDX-License-Identifier: LGPL-3.0+
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

template <class EigenSolver, class FT>
struct Get_eigen_matrix
{
  typedef Eigen_sparse_matrix<FT>             type;
};

template <class FT, class EigenMatrix>
struct Get_eigen_matrix< ::Eigen::ConjugateGradient<EigenMatrix>, FT>
{
  typedef Eigen_sparse_symmetric_matrix<FT>   type;
};

template <class FT, class EigenMatrix>
struct Get_eigen_matrix< ::Eigen::SimplicialCholesky<EigenMatrix>, FT>
{
  typedef Eigen_sparse_symmetric_matrix<FT>   type;
};

#if EIGEN_VERSION_AT_LEAST(3, 1, 91)
template <class FT, class EigenMatrix, class EigenOrdering>
struct Get_eigen_matrix< ::Eigen::SparseLU<EigenMatrix, EigenOrdering >, FT>
{
  typedef Eigen_sparse_matrix<FT>             type;
};
#endif
} //internal

/*!
\ingroup PkgSolver

The class `Eigen_solver_traits` provides an interface to the sparse solvers of \ref thirdpartyEigen "Eigen".
\ref thirdpartyEigen "Eigen" version 3.1 (or later) must be available on the system.

\cgalModels `SparseLinearAlgebraWithFactorTraits_d` and `NormalEquationSparseLinearAlgebraTraits_d`

\tparam EigenSolverT A sparse solver of \ref thirdpartyEigen "Eigen". The default solver is the iterative bi-congugate gradient stabilized solver  `Eigen::BiCGSTAB` for `double`.

\sa `CGAL::Eigen_sparse_matrix<T>`
\sa `CGAL::Eigen_sparse_symmetric_matrix<T>`
\sa `CGAL::Eigen_vector<T>`
\sa http://eigen.tuxfamily.org

Example
-------------- 

The instantiation of this class assumes an \ref thirdpartyEigen "Eigen" sparse solver is provided. Here are few examples:

\code{.cpp}

typedef CGAL::Eigen_sparse_matrix<double>::EigenType EigenMatrix;

//iterative general solver
typedef CGAL::Eigen_solver_traits< Eigen::BiCGSTAB<EigenMatrix> > Iterative_general_solver;

//iterative symmetric solver
typedef CGAL::Eigen_solver_traits< Eigen::ConjugateGradient<EigenMatrix> > Iterative_symmetric_solver;

//direct symmetric solver
typedef CGAL::Eigen_solver_traits< Eigen::SimplicialCholesky<EigenMatrix> > Direct_symmetric_solver;

\endcode
*/
template<class EigenSolverT = Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType> >
class Eigen_solver_traits
{
  typedef typename EigenSolverT::Scalar                               Scalar;

  // Public types
public:
  /// \name Types
  /// @{

  typedef EigenSolverT                                                Solver;
  typedef Scalar                                                      NT;
  typedef CGAL::Eigen_vector<NT>                                      Vector;

  /// If `T` is `Eigen::ConjugateGradient<M>` or `Eigen::SimplicialCholesky<M>`,
  /// `Matrix` is `CGAL::Eigen_sparse_symmetric_matrix<T>`, and `CGAL::Eigen_sparse_matrix<T>` otherwise.
#ifdef DOXYGEN_RUNNING
  typedef unspecified_type                                            Matrix;
#else
  typedef typename internal::Get_eigen_matrix<EigenSolverT,NT>::type  Matrix;
#endif
  /// @}

  // Public operations
public:
  /// Constructor
  Eigen_solver_traits() : m_mat(NULL), m_solver_sptr(new EigenSolverT) { }

  /// \name Operations
  /// @{

  /// Returns a reference to the internal \ref thirdpartyEigen "Eigen" solver.
  /// This function can be used for example to set specific parameters of the solver.
  EigenSolverT& solver() { return *m_solver_sptr; }

  /// @}

  /// Solve the sparse linear system \f$ A \times X = B \f$.
  /// Return `true` on success. The solution is then \f$ (1/D) \times X \f$.
  ///
  /// \pre A.row_dimension() == B.dimension().
  /// \pre A.column_dimension() == X.dimension().
  bool linear_solver(const Matrix& A, const Vector& B, Vector& X, NT& D)
  {
    D = 1; // Eigen does not support homogeneous coordinates

    m_solver_sptr->compute(A.eigen_object());

    if(m_solver_sptr->info() != Eigen::Success)
      return false;

    X = m_solver_sptr->solve(B);

    return m_solver_sptr->info() == Eigen::Success;
  }

  /// Factorize the sparse matrix \f$ A \f$.
  /// This factorization is used in `linear_solver()`
  /// to solve the system for different right-hand side vectors.
  /// See `linear_solver()` for the description of \f$ D \f$.
  /// \return `true` if the factorization is successful and `false` otherwise.
  bool factor(const Matrix& A, NT& D)
  {
    D = 1;

    m_mat = &A.eigen_object();
    solver().compute(*m_mat);
    return solver().info() == Eigen::Success;
  }

  /// Solve the sparse linear system \f$ A \times X = B\f$, with \f$ A \f$ being the matrix
  /// provided in `factor()`.
  /// \return `true` if the solver is successful and `false` otherwise.
  bool linear_solver(const Vector& B, Vector& X)
  {
    CGAL_precondition(m_mat != NULL); // factor should have been called first
    X = solver().solve(B);
    return solver().info() == Eigen::Success;
  }

  /// Factorize the sparse matrix \f$ A^t \times A\f$, where \f$ A^t \f$ is the
  /// transpose matrix of \f$ A \f$.
  /// This factorization is used in `normal_equation_solver()` to solve the system
  /// for different right-hand side vectors.
  /// \return `true` if the factorization is successful and `false` otherwise.
  bool normal_equation_factor(const Matrix& A)
  {
    typename Matrix::EigenType At = A.eigen_object().transpose();
    m_mat = &A.eigen_object();
    solver().compute(At * A.eigen_object());
    return solver().info() == Eigen::Success;
  }

  /// Solve the sparse linear system \f$ A^t \times A \times X = A^t \times B \f$,
  /// with \f$ A \f$ being the matrix provided in `#normal_equation_factor()`,
  /// and \f$ A^t \f$ its transpose matrix.
  /// \return `true` if the solver is successful and `false` otherwise.
  bool normal_equation_solver(const Vector& B, Vector& X)
  {
    CGAL_precondition(m_mat != NULL); // non_symmetric_factor should have been called first
    typename Vector::EigenType AtB = m_mat->transpose() * B.eigen_object();
    X = solver().solve(AtB);
    return solver().info() == Eigen::Success;
  }

  /// Equivalent to a call to \link normal_equation_factor() `normal_equation_factor(A)` \endlink
  /// followed by a call to \link normal_equation_solver `normal_equation_solver(B, X)` \endlink .
  bool normal_equation_solver(const Matrix& A, const Vector& B, Vector& X)
  {
    if (!normal_equation_factor(A))
      return false;
    return normal_equation_solver(B, X);
  }

protected:
  const typename Matrix::EigenType* m_mat;
  boost::shared_ptr<EigenSolverT> m_solver_sptr;
};

// Specialization of the solver for BiCGSTAB as for surface parameterization,
// the intializer should be a vector of one's (this was the case in 3.1-alpha
// but not in the official 3.1).
template<>
class Eigen_solver_traits<Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType> >
{
  typedef Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType>       EigenSolverT;
  typedef EigenSolverT::Scalar                                          Scalar;

  // Public types
public:
  typedef EigenSolverT                                                  Solver;
  typedef Scalar                                                        NT;
  typedef internal::Get_eigen_matrix<EigenSolverT,NT>::type             Matrix;
  typedef Eigen_vector<Scalar>                                          Vector;

  // Public operations
public:
  /// Constructor
  Eigen_solver_traits(): m_solver_sptr(new EigenSolverT) { }

  /// Returns a reference to the internal \ref thirdpartyEigen "Eigen" solver.
  /// This function can be used for example to set specific parameters of the solver.
  EigenSolverT& solver() { return *m_solver_sptr; }

  /// Solve the sparse linear system \f$ A \times X = B \f$.
  /// Return `true` on success. The solution is then \f$ (1/D) \times X \f$.
  ///
  /// \pre A.row_dimension() == B.dimension().
  /// \pre A.column_dimension() == X.dimension().
  bool linear_solver(const Matrix& A, const Vector& B, Vector& X, NT& D)
  {
    D = 1; // Eigen does not support homogeneous coordinates

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

} // namespace CGAL

#endif // CGAL_EIGEN_SOLVER_TRAITS_H
