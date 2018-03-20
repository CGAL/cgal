// Copyright (c) 2015 GeometryFactory (France).
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
//
// Author(s)     : Ilker O. Yaz

#ifndef CGAL_POLYGON_MESH_PROCESSING_FAIR_H
#define CGAL_POLYGON_MESH_PROCESSING_FAIR_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Polygon_mesh_processing/internal/fair_impl.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>  // for sparse linear system solver
#endif

#include <boost/type_traits/is_same.hpp>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {
  // use non-default weight calculator and non-default solver
  // WeightCalculator a model of `FairWeightCalculator`, can be omitted to use default Cotangent weights
  // weight_calculator a function object to calculate weights, defaults to Cotangent weights and can be omitted
  template<typename SparseLinearSolver,
           typename WeightCalculator,
           typename TriangleMesh,
           typename VertexRange,
           typename VertexPointMap>
  bool fair(TriangleMesh& tmesh,
    const VertexRange& vertices,
    SparseLinearSolver solver,
    WeightCalculator weight_calculator,
    unsigned int continuity,
    VertexPointMap vpmap)
  {
    CGAL::Polygon_mesh_processing::internal::Fair_Polyhedron_3
       <TriangleMesh, SparseLinearSolver, WeightCalculator, VertexPointMap>
       fair_functor(tmesh, vpmap, weight_calculator);
    return fair_functor.fair(vertices, solver, continuity);
  }

} //end namespace internal

  /*!
  \ingroup PMP_meshing_grp
  @brief fairs a region on a triangle mesh.
  The points of the selected vertices are
  relocated to yield an as-smooth-as-possible surface patch,
  based on solving a linear bi-Laplacian system with boundary constraints,
  described in \cgalCite{Botsch2008OnLinearVariational}.
  The optional parameter `fairing_continuity` gives the ability to control the tangential
  continuity C<sup> n</sup> of the output mesh.

  The region described by `vertices` might contain multiple disconnected components.
  Note that the mesh connectivity is not altered in any way,
  only vertex locations get updated.

  Fairing might fail if fixed vertices, which are used as boundary conditions,
  do not suffice to solve constructed linear system.

  Note that if the vertex range to which fairing is applied contains all the vertices of the triangle mesh,
  fairing does not fail, but the mesh gets shrinked to `CGAL::ORIGIN`.

  @tparam TriangleMesh a model of `FaceGraph` and `MutableFaceGraph`
  @tparam VertexRange a range of vertex descriptors of `TriangleMesh`, model of `Range`.
          Its iterator type is `InputIterator`.
  @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"

  @param tmesh the triangle mesh with patches to be faired
  @param vertices the vertices of the patches to be faired (the positions of only those vertices will be changed)
  @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `tmesh`.
        If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
    \cgalParamBegin{fairing_continuity} tangential continuity of the output surface patch. The larger `fairing_continuity` gets, the more fixed vertices are required \cgalParamEnd
    \cgalParamBegin{sparse_linear_solver} an instance of the sparse linear solver used for fairing \cgalParamEnd
  \cgalNamedParamsEnd
  
  @return `true` if fairing is successful, otherwise no vertices are relocated

  @pre `is_triangle_mesh(tmesh)`

  @warning This function involves linear algebra, that is computed using a non-exact floating-point arithmetic.

  @todo accuracy of solvers are not good, for example when there is no boundary condition pre_factor should fail, but it does not.
  */
  template<typename TriangleMesh,
           typename VertexRange,
           typename NamedParameters>
  bool fair(TriangleMesh& tmesh,
            const VertexRange& vertices,
            const NamedParameters& np)
  {
    using boost::get_param;
    using boost::choose_param;

    CGAL_precondition(is_triangle_mesh(tmesh));

#if defined(CGAL_EIGEN3_ENABLED)
  #if EIGEN_VERSION_AT_LEAST(3,2,0)
    typedef CGAL::Eigen_solver_traits<Eigen::SparseLU<
      CGAL::Eigen_sparse_matrix<double>::EigenType, Eigen::COLAMDOrdering<int> >  >
      Default_solver;
  #else
    typedef bool Default_solver;//compilation should crash
      //if no solver is provided and Eigen version < 3.2
  #endif
#else
    typedef bool Default_solver;//compilation should crash
      //if no solver is provided and Eigen version < 3.2
#endif

#if defined(CGAL_EIGEN3_ENABLED)
    CGAL_static_assertion_msg(
      (!boost::is_same<typename GetSolver<NamedParameters, Default_solver>::type, bool>::value) || EIGEN_VERSION_AT_LEAST(3, 2, 0),
      "The function `fair` requires Eigen3 version 3.2 or later.");
#else
    CGAL_static_assertion_msg(
      (!boost::is_same<typename GetSolver<NamedParameters, Default_solver>::type, bool>::value),
      "The function `fair` requires Eigen3 version 3.2 or later.");
#endif

    typedef typename GetVertexPointMap < TriangleMesh, NamedParameters>::type VPMap;
    typedef CGAL::internal::Cotangent_weight_with_voronoi_area_fairing<TriangleMesh, VPMap>
      Default_Weight_calculator;

    VPMap vpmap_ = choose_param(get_param(np, internal_np::vertex_point),
                                get_property_map(vertex_point, tmesh));

    return internal::fair(tmesh, vertices,
      choose_param(get_param(np, internal_np::sparse_linear_solver), Default_solver()),
      choose_param(get_param(np, internal_np::weight_calculator), Default_Weight_calculator(tmesh, vpmap_)),
      choose_param(get_param(np, internal_np::fairing_continuity), 1),
      vpmap_
      );
  }

  template<typename TriangleMesh, typename VertexRange>
  bool fair(TriangleMesh& tmesh, const VertexRange& vertices)
  {
    return fair(tmesh, vertices,
      CGAL::Polygon_mesh_processing::parameters::all_default());
  }

} //end namespace Polygon_mesh_processing

} //end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_POLYGON_MESH_PROCESSING_FAIR_H
