// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  @param tmesh the triangle mesh with patches to be faired
  @param vertices the vertices of the patches to be faired (the positions of only those vertices will be changed)
  @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      should be available for the vertices of `tmesh`.}
    \cgalParamNEnd

    \cgalParamNBegin{fairing_continuity}
      \cgalParamDescription{A value controling the tangential continuity of the output surface patch.
                            The possible values are 0, 1 and 2, refering to the  C<sup>0</sup>, C<sup>1</sup>
                            and C<sup>2</sup> continuity.}
      \cgalParamType{unsigned int}
      \cgalParamDefault{`1`}
      \cgalParamExtra{The larger `fairing_continuity` gets, the more fixed vertices are required.}
    \cgalParamNEnd

    \cgalParamNBegin{sparse_linear_solver}
      \cgalParamDescription{an instance of the sparse linear solver used for fairing}
      \cgalParamType{a class model of `SparseLinearAlgebraWithFactorTraits_d`}
      \cgalParamDefault{If \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available and
                        `CGAL_EIGEN3_ENABLED` is defined, then the following overload of `Eigen_solver_traits`
                        is provided as default value:\n
                        `CGAL::Eigen_solver_traits<Eigen::SparseLU<CGAL::Eigen_sparse_matrix<double>::%EigenType, Eigen::COLAMDOrdering<int> > >`}
    \cgalParamNEnd
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
    using parameters::get_parameter;
    using parameters::choose_parameter;

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

    // Cotangent_weight_with_voronoi_area_fairing has been changed to the version:
    // Cotangent_weight_with_voronoi_area_fairing_secure to avoid imprecisions from
    // the issue #4706 - https://github.com/CGAL/cgal/issues/4706.
    typedef CGAL::internal::Cotangent_weight_with_voronoi_area_fairing_secure<TriangleMesh, VPMap>
      Default_Weight_calculator;

    VPMap vpmap_ = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                    get_property_map(vertex_point, tmesh));

    return internal::fair(tmesh, vertices,
      choose_parameter<Default_solver>(get_parameter(np, internal_np::sparse_linear_solver)),
      choose_parameter(get_parameter(np, internal_np::weight_calculator), Default_Weight_calculator(tmesh, vpmap_)),
      choose_parameter(get_parameter(np, internal_np::fairing_continuity), 1),
      vpmap_);
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
