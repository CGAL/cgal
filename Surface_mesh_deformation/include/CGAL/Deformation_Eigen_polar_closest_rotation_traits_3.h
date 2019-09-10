// Copyright (c) 2013  INRIA Bordeaux Sud-Ouest (France), All rights reserved.
// Copyright (c) 2013 GeometryFactory
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
// Author(s)     : Gael Guennebaud Ilker O. Yaz


#ifndef CGAL_DEFORMATION_EIGEN_POLAR_CLOSEST_ROTATION_TRAITS_3_H
#define CGAL_DEFORMATION_EIGEN_POLAR_CLOSEST_ROTATION_TRAITS_3_H

#include <CGAL/Deformation_Eigen_closest_rotation_traits_3.h>
#include <CGAL/FPU_extension.h>
#include <CGAL/Profile_counter.h>

namespace CGAL {
  /// \ingroup PkgSurfaceMeshDeformation
  /// A class to compute the closest rotation in Frobenius norm to a 3x3 Matrix using the \link thirdpartyEigen `Eigen` library \endlink.
  /// The internal computation relies on a hybrid system using the solvers `Eigen::SelfAdjointEigenSolver<>`
  /// and `Eigen::JacobiSVD<>` (polar decomposition).
  ///
  /// \cgalModels `DeformationClosestRotationTraits_3`
  class Deformation_Eigen_polar_closest_rotation_traits_3 :
    public Deformation_Eigen_closest_rotation_traits_3{
  public:

    /// \cond SKIP_FROM_MANUAL

    /// Computes closest rotation to `m` and places it into `R`
    void compute_close_rotation(const Matrix& m, Matrix& R)
    {
      CGAL_PROFILER(" times closest rotation is computed");
      bool solved = polar_eigen(m, R);

      if(!solved) {
        CGAL_PROFILER(" times polar_eigen failed and SVD is called");
        Deformation_Eigen_closest_rotation_traits_3::compute_close_rotation(m, R);
      }
    }

  private:
    // polar decomposition using Eigen, 5 times faster than SVD
    bool polar_eigen(const Matrix& A, Matrix& R)
    {
      if(A.determinant() < 0)
      { return false; }

      typedef Matrix::Scalar Scalar;

      const Scalar th = std::sqrt(Eigen::NumTraits<Scalar>::dummy_precision());

      Eigen::SelfAdjointEigenSolver<Matrix> eig;
      CGAL::feclearexcept(FE_UNDERFLOW);
      eig.computeDirect(A.transpose()*A);
      if(CGAL::fetestexcept(FE_UNDERFLOW) || eig.eigenvalues()(0)/eig.eigenvalues()(2)<th)
      { return false; }

      Vector S = eig.eigenvalues().cwiseSqrt();
      R = A  * eig.eigenvectors() * S.asDiagonal().inverse()
        * eig.eigenvectors().transpose();

      if(std::abs(R.squaredNorm()-3.) > th || R.determinant() < 0)
      { return false; }

      R.transposeInPlace(); // the optimal rotation matrix should be transpose of decomposition result
      return true;
    }
    /// \endcond

  };
}//namespace CGAL
#endif // CGAL_DEFORMATION_EIGEN_POLAR_CLOSEST_ROTATION_TRAITS_3_H
