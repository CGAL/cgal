
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

namespace CGAL{

/// \ingroup PkgVMASRef
/// @brief extracts medial skeleton for the triangle mesh `tmesh`.
/// This function uses the class CGAL::Variational_medial_axis_sampling with the default parameters.
/// @pre `tmesh` is a triangle mesh without borders
/// @pre The specialization `boost::property_map<TriangleMesh, CGAL::vertex_point_t>::%type` and `get(vertex_point, tmesh)` are defined.
/// @pre The value type of `boost::property_map<TriangleMesh, CGAL::vertex_point_t>::%type` is a point type from a \cgal Kernel.
///
/// @tparam TriangleMesh
///         a model of `FaceListGraph`
///
/// @param tmesh
///        input mesh

/// @pre `tmesh` is a triangulated surface mesh without borders 
template <class TriangleMesh, class NamedParameters = parameters::Default_named_parameters>
typename CGAL::Medial_Skeleton<TriangleMesh>
extract_variational_medial_skeleton(const TriangleMesh& tmesh, const NamedParameters& np = parameters::default_values())
{
  using VMAS =  CGAL::Variational_medial_axis<TriangleMesh>;
  VMAS vmas(tmesh);
  vmas.compute_variational_medial_axis_sampling(np);
  return vmas.export_skeleton();
}

}// end of namespace CGAL

#endif //CGAL_EXTRACT_VARIATIONAL_MEDIAL_SKELETON_H
