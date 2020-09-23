// Copyright (c) 2014  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_EXTRACT_MEAN_CURVATURE_FLOW_SKELETON_H
#define CGAL_EXTRACT_MEAN_CURVATURE_FLOW_SKELETON_H

#include <CGAL/license/Surface_mesh_skeletonization.h>


#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_solver_traits.h>  // for sparse linear system solver
#endif

#include <CGAL/Mean_curvature_flow_skeletonization.h>

namespace CGAL{

#if defined(DOXYGEN_RUNNING) || defined(CGAL_EIGEN3_ENABLED)
/// \ingroup PkgSurfaceMeshSkeletonizationRef
/// @brief extracts a medially centered curve skeleton for the triangle mesh `tmesh`.
/// This function uses the class CGAL::Mean_curvature_flow_skeletonization with the default parameters.
/// This function is provided only if \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available and `CGAL_EIGEN3_ENABLED` is defined.
/// @pre `tmesh` is a triangle mesh without borders and having exactly one connected component.
/// @pre The specialization `boost::property_map<TriangleMesh, CGAL::vertex_point_t>::%type` and `get(vertex_point, tmesh)` are defined.
/// @pre The value type of `boost::property_map<TriangleMesh, CGAL::vertex_point_t>::%type` is a point type from a \cgal Kernel.
///
/// @tparam TriangleMesh
///         a model of `FaceListGraph`
///
/// @param tmesh
///        input mesh
/// @param skeleton
///        graph that will contain the skeleton of `tmesh`.
///        For each vertex descriptor `vd` of `skeleton`, the corresponding point
///        and the set of input vertices that contracted to `vd` can be retrieved
///        using `skeleton[vd].point` and `skeleton[vd].vertices` respectively.
///
template <class TriangleMesh>
void extract_mean_curvature_flow_skeleton(const TriangleMesh& tmesh,
                                          typename Mean_curvature_flow_skeletonization<TriangleMesh>::Skeleton& skeleton)
{
  // extract the skeleton
  typedef CGAL::Mean_curvature_flow_skeletonization<TriangleMesh> Mcfskel;
  Mcfskel mcfs(tmesh);
  mcfs(skeleton);
}
#endif

}// end of namespace CGAL

#endif //CGAL_EXTRACT_MEAN_CURVATURE_FLOW_SKELETON_H
