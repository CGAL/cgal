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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Jocelyn Meyron and Quentin Mérigot
//

#ifndef CGAL_VCM_ESTIMATE_EDGES_H
#define CGAL_VCM_ESTIMATE_EDGES_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/vcm_estimate_normals.h>

namespace CGAL {

/// \ingroup PkgPointSetProcessingAlgorithms
/// determines if a point is on a sharp feature edge from a point set
/// for which the Voronoi covariance Measures have been computed.
///
/// The sharpness of the edge, specified by parameter `threshold`, 
/// is used to filtered points according to the external angle around a sharp feature.
///
/// A point is considered to be on a sharp feature if the external angle `alpha` at the edge is such that
/// `alpha >= 2 / sqrt(3) * sqrt(threshold)`.
/// In particular this means that if the input contains sharp features
/// with different external angles, the one with the smallest external angle should be considered,
/// which however would result in selecting more points on sharper regions.
/// More details are provided in \cgalCite{cgal:mog-vbcfe-11}.
///
/// \tparam VCMTraits is a model of `DiagonalizeTraits`. It can be
/// omitted: if Eigen 3 (or greater) is available and
/// `CGAL_EIGEN3_ENABLED` is defined then an overload using
/// `Eigen_diagonalize_traits` is provided. Otherwise, the internal
/// implementation `Diagonalize_traits` is used.
/// \sa CGAL::compute_vcm()`
///
template <class FT, class VCMTraits>
bool
vcm_is_on_feature_edge (cpp11::array<FT,6> &cov,
                        double threshold,
                        VCMTraits)
{
    cpp11::array<double,3> eigenvalues;
    if (!VCMTraits::
          diagonalize_selfadjoint_covariance_matrix(cov, eigenvalues) )
    {
        return false;
    }

    // Compute the ratio
    double r = eigenvalues[1] / (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);
    if (r >= threshold)
        return true;

    return false;
}



template <class FT>
bool
vcm_is_on_feature_edge (cpp11::array<FT,6> &cov,
                        double threshold)
{
  return vcm_is_on_feature_edge(cov, threshold,
				CGAL::Default_diagonalize_traits<double, 3>());

}

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_VCM_ESTIMATE_EDGES_H
