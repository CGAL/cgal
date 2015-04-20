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

#ifndef CGAL_INTERNAL_VCM_VCM_UTILITIES_H
#define CGAL_INTERNAL_VCM_VCM_UTILITIES_H

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace CGAL {

namespace internal {

// Construct the covariance matrix
template <class Covariance>
Eigen::Matrix3f
construct_covariance_matrix (Covariance &cov) {
    Eigen::Matrix3f m;

    m(0,0) = cov[0]; m(0,1) = cov[1]; m(0,2) = cov[2];
    m(1,1) = cov[3]; m(1,2) = cov[4];
    m(2,2) = cov[5];

    m(1, 0) = m(0,1); m(2, 0) = m(0, 2); m(2, 1) = m(1, 2);

    return m;
}

// Diagonalize a selfadjoint matrix
bool
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

// Extract the eigenvector associated to the greatest eigenvalue
template <class Covariance>
bool
extract_greater_eigenvector (Covariance &cov,
                             Eigen::Vector3f &normal) {
    // Construct covariance matrix
    Eigen::Matrix3f m = construct_covariance_matrix(cov);

    // Diagonalizing the matrix
    Eigen::Vector3f eigenvalues;
    Eigen::Matrix3f eigenvectors;
    if (! diagonalize_selfadjoint_matrix(m, eigenvectors, eigenvalues)) {
        return false;
    }

    // Eigenvalues are already sorted by increasing order
    normal = eigenvectors.col(0);

    return true;
}

} // namespace internal

} // namespace CGAL

#endif // CGAL_INTERNAL_VCM_VCM_UTILITIES_H
