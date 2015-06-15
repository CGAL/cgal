// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
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
// Author(s) : Jocelyn Meyron and Quentin Mérigot
//

#ifndef CGAL_EIGEN_VCM_TRAITS_H
#define CGAL_EIGEN_VCM_TRAITS_H

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <CGAL/array.h>

namespace CGAL {

/// A model of the concept `VCMTraits` using \ref thirdpartyEigen.
/// \cgalModels `VCMTraits`
class Eigen_vcm_traits{
  // Construct the covariance matrix
  static Eigen::Matrix3f
  construct_covariance_matrix (const cpp11::array<double,6>& cov) {
    Eigen::Matrix3f m;

    m(0,0) = cov[0]; m(0,1) = cov[1]; m(0,2) = cov[2];
    m(1,1) = cov[3]; m(1,2) = cov[4];
    m(2,2) = cov[5];

    m(1, 0) = m(0,1); m(2, 0) = m(0, 2); m(2, 1) = m(1, 2);

    return m;
  }

  // Diagonalize a selfadjoint matrix
  static bool
  diagonalize_selfadjoint_matrix (Eigen::Matrix3f &m,
                                  Eigen::Matrix3f &eigenvectors,
                                  Eigen::Vector3f &eigenvalues) {
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver(m);

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
    const cpp11::array<double,6>& cov,
    cpp11::array<double, 3>& eigenvalues)
  {
    Eigen::Matrix3f m = construct_covariance_matrix(cov);

    // Diagonalizing the matrix
    Eigen::Vector3f eigenvalues_;
    Eigen::Matrix3f eigenvectors_;
    bool res = diagonalize_selfadjoint_matrix(m, eigenvectors_, eigenvalues_);

    if (res)
    {
      eigenvalues[0]=eigenvalues_[0];
      eigenvalues[1]=eigenvalues_[1];
      eigenvalues[2]=eigenvalues_[2];
    }

    return res;
  }

  // Extract the eigenvector associated to the largest eigenvalue
  static bool
  extract_largest_eigenvector_of_covariance_matrix (
    const cpp11::array<double,6>& cov,
    cpp11::array<double,3> &normal)
  {
      // Construct covariance matrix
      Eigen::Matrix3f m = construct_covariance_matrix(cov);

      // Diagonalizing the matrix
      Eigen::Vector3f eigenvalues;
      Eigen::Matrix3f eigenvectors;
      if (! diagonalize_selfadjoint_matrix(m, eigenvectors, eigenvalues)) {
          return false;
      }

      // Eigenvalues are already sorted by increasing order
      normal[0]=eigenvectors(0,0);
      normal[1]=eigenvectors(1,0);
      normal[2]=eigenvectors(2,0);

      return true;
  }
};

} // namespace CGAL

#endif // CGAL_EIGEN_VCM_TRAITS_H
