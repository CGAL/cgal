// Copyright (c) 2012  INRIA Bordeaux Sud-Ouest (France), All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Gael Guennebaud

#ifndef CGAL_EIGEN_SVD_H
#define CGAL_EIGEN_SVD_H

#include <boost/config.hpp>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4244)
#endif

#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>
#include <Eigen/SVD>

namespace CGAL {

  /*!
\ingroup PkgSolverInterfaceRef

The class `Eigen_svd` provides an algorithm to solve in the least
square sense a linear system with a singular value decomposition using
\ref thirdpartyEigen.

\cgalModels `SvdTraits`

*/
class Eigen_svd
{
public:
  /// \name Types
  /// @{

  typedef double                                          FT;
  typedef Eigen_vector<FT>                                Vector;
  typedef Eigen_matrix<FT>                                Matrix;

  /// @}

  /// Solves the system \f$ MX=B\f$ (in the least square sense if \f$ M\f$ is not
  /// square) using a singular value decomposition.The solution is stored in \f$ B\f$.
  /// \return the condition number of \f$ M\f$
  static FT solve(const Matrix& M, Vector& B)
  {
    Eigen::JacobiSVD<Matrix::EigenType> jacobiSvd(M.eigen_object(),::Eigen::ComputeThinU | ::Eigen::ComputeThinV);
    B.eigen_object()=jacobiSvd.solve(Vector::EigenType(B.eigen_object()));
    return jacobiSvd.singularValues().array().abs().maxCoeff() /
           jacobiSvd.singularValues().array().abs().minCoeff();
  }
};

} // namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_EIGEN_SVD_H
