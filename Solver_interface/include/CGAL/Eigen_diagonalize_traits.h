// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// Author(s) : Jocelyn Meyron and Quentin Mérigot
//

#ifndef CGAL_EIGEN_DIAGONALIZE_TRAITS_H
#define CGAL_EIGEN_DIAGONALIZE_TRAITS_H

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


// If the matrix to diagonalize is of dimension 2x2 or 3x3, Eigen
// provides a faster implementation using a closed-form
// algorithm. However, it offers less precision. See:
// https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html
// This is usually acceptable for CGAL algorithms but one might want
// to use the slower but more accurate version. In that case, just
// uncomment the following line:

//#define DO_NOT_USE_EIGEN_COMPUTEDIRECT_FOR_DIAGONALIZATION

#include <CGAL/array.h>

namespace CGAL {

/// A model of the concept `DiagonalizeTraits` using \ref thirdpartyEigen.
/// \cgalModels `DiagonalizeTraits`

template <typename FT, unsigned int dim = 3>
class Eigen_diagonalize_traits{

public:

  typedef cpp11::array<FT, dim> Vector;
  typedef cpp11::array<FT, dim*dim> Matrix;
  typedef cpp11::array<FT, (dim * (dim+1) / 2)> Covariance_matrix;
  
private:

  typedef Eigen::Matrix<FT, dim, dim> EigenMatrix;
  typedef Eigen::Matrix<FT, dim, 1> EigenVector;
  
  // Construct the covariance matrix
  static EigenMatrix
  construct_covariance_matrix
  (const Covariance_matrix& cov)  {
    EigenMatrix m;

    for (std::size_t i = 0; i < dim; ++ i)
      for (std::size_t j = i; j < dim; ++ j)
	{
	  m(i,j) = static_cast<float>(cov[(dim * i) + j - ((i * (i+1)) / 2)]);
	  if (i != j)
	    m(j,i) = m(i,j);
	}

    return m;
  }

  // Diagonalize a selfadjoint matrix
  static bool
  diagonalize_selfadjoint_matrix (EigenMatrix& m,
				  EigenMatrix& eigenvectors,
                                  EigenVector& eigenvalues) {
      Eigen::SelfAdjointEigenSolver<EigenMatrix> eigensolver;

#ifndef DO_NOT_USE_EIGEN_COMPUTEDIRECT_FOR_DIAGONALIZATION
      if (dim == 2 || dim == 3)
        eigensolver.computeDirect(m);
      else
#endif
        eigensolver.compute(m);

      if (eigensolver.info() != Eigen::Success) {
          return false;
      }

      eigenvalues = eigensolver.eigenvalues();
      eigenvectors = eigensolver.eigenvectors();

      return true;
  }

public:

  static bool
  diagonalize_selfadjoint_covariance_matrix(
    const Covariance_matrix& cov,
    Vector& eigenvalues)
  {
    EigenMatrix m = construct_covariance_matrix(cov);

    // Diagonalizing the matrix
    EigenVector eigenvalues_;
    EigenMatrix eigenvectors_;
    bool res = diagonalize_selfadjoint_matrix(m, eigenvectors_, eigenvalues_);

    if (res)
    {
      for (std::size_t i = 0; i < dim; ++ i)
	eigenvalues[i] = static_cast<FT>(eigenvalues_[i]);
    }

    return res;
  }

  static bool
  diagonalize_selfadjoint_covariance_matrix(
    const Covariance_matrix& cov,
    Vector& eigenvalues,
    Matrix& eigenvectors)
  {
    EigenMatrix m = construct_covariance_matrix(cov);

    // Diagonalizing the matrix
    EigenVector eigenvalues_;
    EigenMatrix eigenvectors_;
    bool res = diagonalize_selfadjoint_matrix(m, eigenvectors_, eigenvalues_);

    if (res)
    {
      for (std::size_t i = 0; i < dim; ++ i)
	{
	  eigenvalues[i] = static_cast<FT>(eigenvalues_[i]);

	  for (std::size_t j = 0; j < dim; ++ j)
	    eigenvectors[dim*i + j]=static_cast<FT>(eigenvectors_(j,i));
	}
    }

    return res;
  }

  // Extract the eigenvector associated to the largest eigenvalue
  static bool
  extract_largest_eigenvector_of_covariance_matrix (
    const Covariance_matrix& cov,
    Vector& normal)
  {
      // Construct covariance matrix
      EigenMatrix m = construct_covariance_matrix(cov);

      // Diagonalizing the matrix
      EigenVector eigenvalues;
      EigenMatrix eigenvectors;
      if (! diagonalize_selfadjoint_matrix(m, eigenvectors, eigenvalues)) {
          return false;
      }

      // Eigenvalues are sorted by increasing order
      for (unsigned int i = 0; i < dim; ++ i)
	normal[i] = static_cast<FT> (eigenvectors(i,dim-1));

      return true;
  }
};

} // namespace CGAL

#endif // CGAL_EIGEN_DIAGONALIZE_TRAITS_H
