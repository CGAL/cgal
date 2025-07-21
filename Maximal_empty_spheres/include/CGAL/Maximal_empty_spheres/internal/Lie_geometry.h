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

#ifndef CGAL_MAXIMAL_EMPTY_SPHERES_INTERNAL_LIE_GEOMETRY_H
#define CGAL_MAXIMAL_EMPTY_SPHERES_INTERNAL_LIE_GEOMETRY_H

#include <CGAL/license/Maximal_empty_spheres.h>


#include <Eigen/Core>

namespace CGAL {
namespace Maximal_empty_spheres {
namespace internal {

inline void lie_ip_matrix(const int D, Eigen::MatrixXd &H){
    H.resize(D+3, D+3);
    H.setZero();
    H.block(0,0,D,D).setIdentity();
    H(D  ,D+1) = -1.;
    H(D+1,D  ) = -1.;
    H(D+2,D+2) = -1.;
}

inline void spheres_to_lie(const Eigen::MatrixXd &C, Eigen::MatrixXd &S){
    const int N = C.rows();
    const int D = C.cols() - 1;

    S.resize(C.rows(),D+3);
    S.setZero();

    S.block(0,0,C.rows(),D) = C.block(0,0,C.rows(),D);

    S.col(D) = 0.5 * (C.block(0,0,C.rows(),D).array().square().rowwise().sum() - C.col(D).array().square());
    S.col(D+1) = Eigen::VectorXd::Ones(N);
    S.col(D+2) = C.col(D);
}

inline void lie_to_spheres(const Eigen::MatrixXd &S, Eigen::MatrixXd &C){
    // assumes x_lie to lie on the lie quadric and have finite radius
    const int N = S.rows();
    const int D = S.cols() - 3;

    C.resize(N,D+1);
    C.block(0,0,N,D) = S.block(0,0,N,D).array().colwise() / S.col(D+1).array();
    C.col(D) = S.col(D+2).array() / S.col(D+1).array();
}

inline void line_quadric_intersection(const Eigen::VectorXd &x, const Eigen::VectorXd &y, const Eigen::MatrixXd &H, double &lambda_1, double &lambda_2){
    // intersect line (1-lambda) x + lambda y with the lie quadric x_^T H x_ = 0

    double a =     (y-x).transpose() * H * (y-x);
    double b = 2 * (y-x).transpose() * H * x;
    double c =         x.transpose() * H * x;

    double delta = b*b-4*a*c;
    lambda_1 = (-b-sqrt(delta)) / (2*a);
    lambda_2 = (-b+sqrt(delta)) / (2*a);
}

} // namespace internal
} // namespace Maximal_empty_spheres
} // namespace CGAL
#endif // CGAL_MAXIMAL_EMPTY_SPHERES_INTERNAL_LIE_GEOMETRY_H
