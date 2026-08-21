// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Qijia Huang

#ifndef CGAL_EXTRACT_VARIATIONAL_MEDIAL_SKELETON_H
#define CGAL_EXTRACT_VARIATIONAL_MEDIAL_SKELETON_H

#include <CGAL/license/Variational_medial_axis.h>

#include <CGAL/Variational_medial_axis_sampling.h>

namespace CGAL {

/// \ingroup PkgVMASRef
/// @brief extracts a medial skeleton for the triangle mesh `tmesh`.
/// This function uses the class CGAL::Variational_medial_axis_sampling with the default parameters.
/// @pre `tmesh` is a triangle mesh without borders
/// @pre The specialization `boost::property_map<TriangleMesh, CGAL::vertex_point_t>::%const_type` and `get(CGAL::vertex_point, tmesh)` are defined.
///
/// @pre The value type of `boost::property_map<TriangleMesh, CGAL::vertex_point_t>::%const_type` is a point type from a \cgal
/// Kernel.
///
/// @tparam TriangleMesh
///         a model of `FaceListGraph`
/// @tparam NamedParameters
///         a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// @param tmesh
///        input mesh
/// @param np
///        an optional sequence of \ref bgl_namedparameters "Named Parameters", listed below:
///
/// \cgalNamedParamsBegin
///  \cgalParamNBegin{number_of_spheres}
///    \cgalParamDescription{The desired number of medial spheres in the resulting skeleton.
///                          Note that this number might not be achieved if the maximum number of iterations is reached or if
///                          the internal energy minimisation converged before. If not provided, the number of spheres will be
///                          the one at the end of the energy minimisation or at the end of the iterations.}
///    \cgalParamType{unsigned int}
///    \cgalParamExtra{This number should generally not exceed 300, as the method is designed to produce coarse skeletons.}
///  \cgalParamNEnd
///   \cgalParamNBegin{neighbor_radius}
///     \cgalParamDescription{This parameter is used to sample the input triangle mesh for the optimization process.
///                           Sample points will be chosen such that it should be possible to cover the surface area
///                           of the triangle mesh using circle of that radius with centered at the sample points.}
///     \cgalParamType{unsigned int}
///     \cgalParamDefault{A radius corresponding to 20,000 sample points}
///  \cgalParamNEnd
///  \cgalParamNBegin{number_of_iterations}
///     \cgalParamDescription{The maximum number of iterations for the optimization process.}
///     \cgalParamType{int}
///     \cgalParamDefault{1000}
///     \cgalParamExtra{This parameter must be strictly positive; setting it to zero may prevent correct skeleton
/// connectivity construction.}
///   \cgalParamNEnd
///  \cgalParamNBegin{lambda}
///    \cgalParamDescription{A weight balancing the two energy terms (SQEM and Euclidean). Smaller values tend to
/// produce skeletons that follow local features more closely.}
///    \cgalParamType{FT}
///    \cgalParamDefault{FT(0.2)}
///    \cgalParamExtra{The range of this parameter is (0,1].}
///  \cgalParamNEnd
///   \cgalParamNBegin{random_seed}
///     \cgalParamDescription{The random seed to sample points on the triangle mesh surface.}
///     \cgalParamType{unsigned int}
///     \cgalParamExtra{Fix the random seed so that the result can be reproduced}
///   \cgalParamNEnd
///  \cgalParamNBegin{concurrency_tag}
///    \cgalParamDescription{Tag indicating whether the algorithm should run sequentially or in parallel.}
///    \cgalParamType{Concurrency tag type}
///    \cgalParamDefault{`CGAL::Sequential_tag`}
///    \cgalParamExtra{Use `CGAL::Parallel_tag` for parallel execution (requires TBB).}
///  \cgalParamNEnd
/// \cgalNamedParamsEnd
///
///
template <class TriangleMesh, class NamedParameters = parameters::Default_named_parameters>
typename CGAL::Medial_skeleton<TriangleMesh>
extract_variational_medial_skeleton(const TriangleMesh& tmesh,
                                    const NamedParameters& np = parameters::default_values()) {
  typedef typename internal_np::Lookup_named_param_def<internal_np::concurrency_tag_t, NamedParameters,
                                                       Sequential_tag>::type Concurrency_tag;
  using VMAS = CGAL::Variational_medial_axis_sampling<TriangleMesh, Concurrency_tag>;
  VMAS vmas(tmesh, np);
  vmas.sample(np);
  return vmas.skeleton();
}
} // end of namespace CGAL

#endif // CGAL_EXTRACT_VARIATIONAL_MEDIAL_SKELETON_H
