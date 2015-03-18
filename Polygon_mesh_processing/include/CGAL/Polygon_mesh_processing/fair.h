#ifndef CGAL_POLYGON_MESH_PROCESSING_FAIR_H
#define CGAL_POLYGON_MESH_PROCESSING_FAIR_H

#include <CGAL/Polygon_mesh_processing/internal/fair_impl.h>

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
    unsigned int continuity = 1)
  {
    CGAL::Polygon_mesh_processing::internal::Fair_Polyhedron_3<PolygonMesh,
      SparseLinearSolver, WeightCalculator> fair_functor(pmesh, weight_calculator);
    return fair_functor.fair(vertices, solver, continuity);
  }

#if defined(CGAL_EIGEN3_ENABLED) && EIGEN_VERSION_AT_LEAST(3,2,0)
  //use default SparseLinearSolver (used in triangulate_hole.h)
  template<class WeightCalculator, class PolygonMesh, class VertexRange>
  bool fair(PolygonMesh& pmesh,
    VertexRange vertices,
    CGAL::Default,
    WeightCalculator weight_calculator,
    unsigned int continuity = 1)
  {
    typedef   CGAL::Eigen_solver_traits<
                Eigen::SparseLU< CGAL::Eigen_sparse_matrix<double>::EigenType,
                                 Eigen::COLAMDOrdering<int> >  >
                                    Sparse_linear_solver;
    return internal::fair(pmesh, vertices, Sparse_linear_solver(), weight_calculator, continuity);
  }
#endif
} //end namespace internal

  /*!
  \ingroup PkgPolygonMeshProcessing
  @brief Fairs a region on a polygon mesh.
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
  template<class SparseLinearSolver, class PolygonMesh, class VertexRange>
  bool fair(PolygonMesh& pmesh,
    VertexRange vertices,
    SparseLinearSolver solver
#ifdef DOXYGEN_RUNNING
    = CGAL::Default()
#endif
    , unsigned int continuity = 1)
  {
    typedef CGAL::internal::Cotangent_weight_with_voronoi_area_fairing<PolygonMesh> Weight_calculator;
    return internal::fair<SparseLinearSolver, Weight_calculator, PolygonMesh, VertexRange>
      (pmesh, vertices, solver, Weight_calculator(pmesh), continuity);
  }

#if defined(CGAL_EIGEN3_ENABLED) && EIGEN_VERSION_AT_LEAST(3,2,0)
  // use default SparseLinearSolver and continuity

  template<class PolygonMesh, class VertexRange>
  bool fair(PolygonMesh& pmesh,
    VertexRange vertices,
    CGAL::Default,
    unsigned int continuity = 1)
  {
    typedef   CGAL::Eigen_solver_traits<
                Eigen::SparseLU< CGAL::Eigen_sparse_matrix<double>::EigenType,
                                 Eigen::COLAMDOrdering<int> >  >
                                    Sparse_linear_solver;
    return fair(pmesh, vertices, Sparse_linear_solver(), continuity);
  }

  template<class PolygonMesh, class VertexRange>
  bool fair(PolygonMesh& pmesh,
            VertexRange vertices)
  {
    return fair(pmesh, vertices, CGAL::Default());
  }
#endif

} //end namespace Polygon_mesh_processing

} //end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_FAIR_H
