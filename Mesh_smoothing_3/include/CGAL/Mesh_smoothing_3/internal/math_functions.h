#pragma once
#include <Eigen/Eigen>

namespace Mesh_smoothing_3_internal {

namespace Math_functions {
    inline Eigen::Matrix<double,1,3> sub_line_vector(Eigen::VectorXd const &x, unsigned i) {
        i *= 3;
        return {x(i), x(i+1), x(i+2)};
    }

    inline Eigen::Matrix<double,3,1> sub_col_vector(Eigen::VectorXd const &x, unsigned i) {
        i *= 3;
        return {x(i), x(i+1), x(i+2)};
    }


    inline double chi(double eps, double det) {
        double const eps2 = eps * eps;
        return det > 0 ? // for numerical stability
                (det + std::sqrt(eps2 + det * det)) * .5 :
                .5 * eps2 / (std::sqrt(eps2 + det * det) - det);
    }

    inline double chi_deriv(double eps, double det) {
        return .5 + det / (2. * std::sqrt(eps * eps + det * det));
    }

    inline std::array<Eigen::Vector3d, 4> transform_coordinates_to_gradient_base(std::array<Eigen::Vector3d, 4> const &vertices_coordinates) {
        Eigen::Matrix3d M;
        M.row(0) = vertices_coordinates[1] - vertices_coordinates[0];
        M.row(1) = vertices_coordinates[2] - vertices_coordinates[0];
        M.row(2) = vertices_coordinates[3] - vertices_coordinates[0];

        Eigen::Matrix3d invM = M.inverse();
        return { -invM.col(0) - invM.col(1) - invM.col(2), invM.col(0), invM.col(1), invM.col(2) };
    }

    // dual basis; i.e. d(detJ)/dJ
    inline Eigen::Matrix3d dual_basis(Eigen::Matrix3d const &J) {
        Eigen::Matrix3d K;
        K.col(0) = J.col(1).cross(J.col(2));
        K.col(1) = J.col(2).cross(J.col(0));
        K.col(2) = J.col(0).cross(J.col(1));
        return K;
    }

}

}

