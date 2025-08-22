// Copyright (c) 2025  TU Berlin
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Max Kohlbrenner

#ifndef CGAL_MAXIMAL_EMPTY_SPHERES_INTERNAL_ROTATION_FROM_TO_H
#define CGAL_MAXIMAL_EMPTY_SPHERES_INTERNAL_ROTATION_FROM_TO_H

// #include <CGAL/license/Maximal_empty_spheres.h>

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SVD>

namespace CGAL {
namespace Maximal_empty_spheres {
namespace internal {

inline Eigen::MatrixXd rotation_from_to(const Eigen::VectorXd &source, const Eigen::VectorXd &target, int debug_level = 0) {
    // get an arbitrary rotation matrix that rotations source to target direction

    if (source.size() != target.size()) {
        std::cout << "(ERROR:) trying to rotate in different dimensions" << std::endl;
    }

    int d = source.size();

    // Rs rotates first unit vector into source
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(d,d);
    S.col(0) = source;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    svd.compute(S, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd Rs = svd.matrixU()*svd.matrixV().transpose();
    if (Rs.determinant() < 1e-10) Rs.col(d-1) *= -1;

    // Rt rotates first unit vectorinto target
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(d,d);
    T.col(0) = target;
    svd.compute(T, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd Rt = svd.matrixU()*svd.matrixV().transpose();
    if (Rt.determinant() < 1e-10) Rt.col(d-1) *= -1;

    // assemble the rotation
    Eigen::MatrixXd R = Rt*Rs.transpose();

    if (debug_level >= 1){
        std::cout << "S: " << std::endl << S << std::endl;
        std::cout << "Rs: " << std::endl << Rs << std::endl;
        std::cout << "T: " << std::endl << T << std::endl;
        std::cout << "Rt: " << std::endl << Rt << std::endl;
        std::cout << "R: " << std::endl << R << std::endl;
    }

    if (debug_level >= 2){
        std::cout << "Rotation det: " << R.determinant() << std::endl;
    }

    return R;
}

} // namespace internal
} // namespace Maximal_empty_spheres
} // namespace CGAL
#endif // CGAL_MAXIMAL_EMPTY_SPHERES_INTERNAL_ROTATION_FROM_TO_H
