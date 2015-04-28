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
// 
//
// Author(s)     : Ilker O. Yaz

#ifndef CGAL_POLYGON_MESH_PROCESSING_FAIR_H
#define CGAL_POLYGON_MESH_PROCESSING_FAIR_H

#include <CGAL/Polygon_mesh_processing/internal/fair_impl.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>  // for sparse linear system solver
#endif


namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {
  // use non-default weight calculator and non-default solver
  // WeightCalculator a model of `FairWeightCalculator`, can be omitted to use default Cotangent weights
  // weight_calculator a function object to calculate weights, defaults to Cotangent weights and can be omitted
  template<typename SparseLinearSolver,
           typename WeightCalculator,
           typename PolygonMesh,
           typename VertexRange,
           typename VertexPointMap>
  bool fair(PolygonMesh& pmesh,
    const VertexRange& vertices,
    SparseLinearSolver solver,
    WeightCalculator weight_calculator,
    unsigned int continuity,
    VertexPointMap vpmap)
  {
    CGAL::Polygon_mesh_processing::internal::Fair_Polyhedron_3
       <PolygonMesh, SparseLinearSolver, WeightCalculator, VertexPointMap>
       fair_functor(pmesh, vpmap, weight_calculator);
    return fair_functor.fair(vertices, solver, continuity);
  }

} //end namespace internal

  /*!
  \ingroup PkgPolygonMeshProcessing
  @brief fairs a region on a polygon mesh.
  The points of the selected vertices are moved to get an as-smooth-as-possible surface patch.
  
  The region described by `vertices` might contain multiple disconnected components.
  Note that the structure is not altered in any way, only positions of the vertices get updated.

  Fairing might fail if fixed vertices, which are used as boundary conditions, do not suffice to solve constructed linear system.

  @tparam PolygonMesh a model of `FaceGraph`
          that has an internal property map for `CGAL::vertex_point_t`
  @tparam VertexRange a range of vertex descriptors of `PolygonMesh`, model of `Range`.
          Its iterator type is `InputIterator`.
  @tparam NamedParameters a sequence of \ref namedparameters

  @param pmesh the polygon mesh with patches to be faired
  @param vertices the vertices of the patches to be faired (the positions of only those vertices will be changed)
  @param np optional sequence of \ref namedparameters among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
    \cgalParamBegin{fairing_continuity} tangential continuity of the output surface patch. The larger `fairing_continuity` gets, the more fixed vertices are required \cgalParamEnd
    \cgalParamBegin{sparse_linear_solver} an instance of the sparse linear solver used for fairing \cgalParamEnd
  \cgalNamedParamsEnd
  
  @return `true` if fairing is successful, otherwise no vertex position is changed

  @todo accuracy of solvers are not good, for example when there is no boundary condition pre_factor should fail, but it does not.
  */
  template<typename PolygonMesh,
           typename VertexRange,
           typename NamedParameters>
  bool fair(PolygonMesh& pmesh,
            const VertexRange& vertices,
            const NamedParameters& np)
  {
    using boost::get_param;
    using boost::choose_param;

    typedef CGAL::internal::Cotangent_weight_with_voronoi_area_fairing<PolygonMesh>
      Default_Weight_calculator;

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

    return internal::fair(pmesh, vertices,
      choose_param(get_param(np, sparse_linear_solver), Default_solver()),
      choose_param(get_param(np, weight_calculator), Default_Weight_calculator(pmesh)),
      choose_param(get_param(np, fairing_continuity), 1),
      choose_param(get_param(np, vertex_point), get(CGAL::vertex_point, pmesh))
      );
  }

  template<typename PolygonMesh, typename VertexRange>
  bool fair(PolygonMesh& pmesh, const VertexRange& vertices)
  {
    return fair(pmesh, vertices,
      CGAL::Polygon_mesh_processing::parameters::all_default());
  }
  
} //end namespace Polygon_mesh_processing

} //end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_FAIR_H
