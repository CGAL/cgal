
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

#include <CGAL/variational_medial_axis_sampling.h>

namespace CGAL {

/// \ingroup PkgVMASRef
/// @brief extracts medial skeleton for the triangle mesh `tmesh`.
/// This function uses the class CGAL::Variational_medial_axis_sampling with the default parameters.
/// @pre `tmesh` is a triangle mesh without borders
/// @pre The specialization `boost::property_map<TriangleMesh, CGAL::vertex_point_t>::%type` and `get(vertex_point, tmesh)` are defined.
///
/// @pre The value type of `boost::property_map<TriangleMesh, CGAL::vertex_point_t>::%type` is a point type from a \cgal
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
///    \cgalParamDescription{The desired number of medial spheres in the resulting skeleton.}
///    \cgalParamType{unsigned int}
///    \cgalParamDefault{100}
///    \cgalParamExtra{This number should generally not exceed 300, as the method is designed to produce coarse
/// skeletons.}
///  \cgalParamNEnd
///  \cgalParamNBegin{max_iteration}
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
///    \cgalParamExtra{This parameter must be strictly positive; setting it to zero may prevent correct skeleton
/// connectivity construction.}
///  \cgalParamNEnd
///  \cgalParamNBegin{concurrency_tag}
///    \cgalParamDescription{Tag indicating whether the algorithm should run sequentially or in parallel.}
///    \cgalParamType{Concurrency tag type}
///    \cgalParamDefault{`CGAL::Sequential_tag`}
///    \cgalParamExtra{Use `CGAL::Parallel_tag` for parallel execution (requires TBB).}
///  \cgalParamNEnd
///  \cgalParamNBegin{acceleration_structure_tag}
///    \cgalParamDescription{Tag indicating the type of acceleration structure to use.}
///    \cgalParamType{Acceleration structure tag type}
///    \cgalParamDefault{`CGAL::KD_tree_tag`}
///    \cgalParamExtra{Use `CGAL::BVH_tag` for a bounding volume hierarchy.}
///  \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// @pre `tmesh` is a triangulated surface mesh without borders
///
template <class TriangleMesh, class NamedParameters = parameters::Default_named_parameters>
typename CGAL::Medial_Skeleton<TriangleMesh>
extract_variational_medial_skeleton(const TriangleMesh& tmesh,
                                    const NamedParameters& np = parameters::default_values()) {
  typedef typename internal_np::Lookup_named_param_def<internal_np::concurrency_tag_t, NamedParameters,
                                                       Sequential_tag>::type Concurrency_tag;
  typedef typename internal_np::Lookup_named_param_def<internal_np::acceleration_structure_t, NamedParameters,
                                                        KD_tree_tag>::type Acceleration_type;
  using VMAS = CGAL::Variational_medial_axis<TriangleMesh, Concurrency_tag, Acceleration_type>;
  VMAS vmas(tmesh, np);
  vmas.compute_variational_medial_axis_sampling(np);
  return vmas.export_skeleton();

}
} // end of namespace CGAL

#endif // CGAL_EXTRACT_VARIATIONAL_MEDIAL_SKELETON_H
