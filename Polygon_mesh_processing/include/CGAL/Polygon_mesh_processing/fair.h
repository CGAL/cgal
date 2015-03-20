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
  template<class SparseLinearSolver, class WeightCalculator, class PolygonMesh, class VertexRange>
  bool fair(PolygonMesh& pmesh,
    VertexRange vertices,
    SparseLinearSolver solver,
    WeightCalculator weight_calculator,
    unsigned int continuity)
  {
    CGAL::Polygon_mesh_processing::internal::Fair_Polyhedron_3<PolygonMesh,
      SparseLinearSolver, WeightCalculator> fair_functor(pmesh, weight_calculator);
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
  The larger `continuity` gets, the more fixed vertices are required.


  @tparam SparseLinearSolver a model of SparseLinearAlgebraTraitsWithFactor_d. If \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available
  and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits` is provided as default parameter:\n
   \code
      CGAL::Eigen_solver_traits<
          Eigen::SparseLU<
             CGAL::Eigen_sparse_matrix<double>::EigenType,
             Eigen::COLAMDOrdering<int> >  >
  \endcode
  @tparam PolygonMesh a model of `FaceGraph`
  @tparam VertexRange a range of vertex descriptors of `PolygonMesh`, model of `SinglePassRange`

  @param pmesh the polygon mesh with patches to be faired
  @param vertices the vertices of the patches to be faired (the positions of only those vertices will be changed)
  @param solver an instance of the sparse linear solver to use.
  @param continuity continuity at the patch boundary (0, 1 and 2 are the possible values)

  @return `true` if fairing is successful, otherwise no vertex position is changed

  \todo SUBMISSION: missing VertexPointMap
  @todo accuracy of solvers are not good, for example when there is no boundary condition pre_factor should fail, but it does not.
  */
  template<typename PolygonMesh,
           typename VertexRange,
           class P, class T, class R>
  bool fair(PolygonMesh& pmesh,
            VertexRange vertices,
            const pmp_bgl_named_params<P, T, R>& p)
  {
    using boost::get_param;
    using boost::choose_param;

    typedef CGAL::internal::Cotangent_weight_with_voronoi_area_fairing<PolygonMesh>
      Default_Weight_calculator;

#if defined(CGAL_EIGEN3_ENABLED) && EIGEN_VERSION_AT_LEAST(3,2,0)
    typedef CGAL::Eigen_solver_traits<Eigen::SparseLU<
      CGAL::Eigen_sparse_matrix<double>::EigenType, Eigen::COLAMDOrdering<int> >  >
      Default_solver;
#else
    typedef bool Default_solver;//compilation should crash
      //if no solver is provided and Eigen version < 3.2
#endif

    return internal::fair(pmesh, vertices,
      choose_param(get_param(p, sparse_linear_solver), Default_solver()),
      choose_param(get_param(p, weight_calculator), Default_Weight_calculator(pmesh)),
      choose_param(get_param(p, fairing_continuity), 1)
      );
  }

  template<typename PolygonMesh, typename VertexRange>
  bool fair(PolygonMesh& pmesh, VertexRange vertices)
  {
    return fair(pmesh, vertices, CGAL::parameters::all_default());
  }
  
} //end namespace Polygon_mesh_processing

} //end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_FAIR_H
