#ifndef CGAL_POLYGON_MESH_PROCESSING_FAIR_H
#define CGAL_POLYGON_MESH_PROCESSING_FAIR_H

#include <CGAL/internal/Meshing_functions/Fair_Polyhedron_3.h>

namespace CGAL {

namespace Polygon_mesh_processing {

  /*!
  \ingroup PkgPolygonMeshProcessing
  @brief Function fairing a region on polygon mesh.
  The region denoted by @a vertex_begin and @a vertex_end might contain multiple disconnected components.
  Note that the structure is not altered in any way, only positions of the vertices get updated.

  Fairing might fail if fixed vertices, which are used as boundary conditions, do not suffice to solve constructed linear system.
  The larger @a continuity gets, the more fixed vertices are required.

  @tparam SparseLinearSolver a model of `SparseLinearAlgebraTraitsWithPreFactor_d`. If \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available
  and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits` is provided as default parameter.
  @tparam PolygonMesh must be a model of `MutableFaceGraph`
  @tparam InputIterator iterator over input vertices

  @param pmesh polygon mesh to be faired
  @param vertex_begin first iterator of the range of vertices
  @param vertex_end past-the-end iterator of the range of vertices
  @param continuity tangential continuity, default to `FAIRING_C_1` and can be omitted

  @return `true` if fairing is successful, otherwise no vertex position is changed

  @todo accuracy of solvers are not good, for example when there is no boundary condition pre_factor should fail, but it does not.
  */
  template<class SparseLinearSolver, class PolygonMesh, class InputIterator>
  bool fair(PolygonMesh& pmesh,
    InputIterator vb,
    InputIterator ve,
    Fairing_continuity continuity = FAIRING_C_1
    )
  {
    typedef CGAL::internal::Cotangent_weight_with_voronoi_area_fairing<PolygonMesh> Weight_calculator;
    return fair<SparseLinearSolver, Weight_calculator, PolygonMesh, InputIterator>
      (pmesh, vb, ve, Weight_calculator(pmesh), continuity);
  }

  // use non-default weight calculator
  // WeightCalculator a model of `FairWeightCalculator`, can be omitted to use default Cotangent weights
  // weight_calculator a function object to calculate weights, defaults to Cotangent weights and can be omitted
  template<class SparseLinearSolver, class WeightCalculator, class PolygonMesh, class InputIterator>
  bool fair(PolygonMesh& pmesh,
    InputIterator vertex_begin,
    InputIterator vertex_end,
    WeightCalculator weight_calculator,
    Fairing_continuity continuity = FAIRING_C_1
    )
  {
    CGAL::Polygon_mesh_processing::internal::Fair_Polyhedron_3<PolygonMesh,
      SparseLinearSolver, WeightCalculator> fair_functor(pmesh, weight_calculator);
    return fair_functor.fair(vertex_begin, vertex_end, continuity);
  }

  //use default SparseLinearSolver
  template<class WeightCalculator, class PolygonMesh, class InputIterator>
  bool fair(PolygonMesh& pmesh,
    InputIterator vb,
    InputIterator ve,
    WeightCalculator weight_calculator,
    Fairing_continuity continuity = FAIRING_C_1
    )
  {
    typedef internal::Fair_default_sparse_linear_solver::Solver Sparse_linear_solver;
    return fair<Sparse_linear_solver, WeightCalculator, PolygonMesh, InputIterator>
      (pmesh, vb, ve, weight_calculator, continuity);
  }

  //use default SparseLinearSolver and WeightCalculator
  template<class PolygonMesh, class InputIterator>
  bool fair(PolygonMesh& pmesh,
    InputIterator vb,
    InputIterator ve,
    Fairing_continuity continuity = FAIRING_C_1
    )
  {
    typedef internal::Fair_default_sparse_linear_solver::Solver Sparse_linear_solver;
    return fair<Sparse_linear_solver, PolygonMesh, InputIterator>(pmesh, vb, ve, continuity);
  }

} //end namespace Polygon_mesh_processing

} //end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_FAIR_H
