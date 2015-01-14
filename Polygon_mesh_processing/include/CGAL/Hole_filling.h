#ifndef CGAL_HOLE_FILLING_H
#define CGAL_HOLE_FILLING_H

// Helper functions which combine triangulate, fair, and refine 

#include <CGAL/internal/Hole_filling/Fair_Polyhedron_3.h>
#include <CGAL/internal/Hole_filling/Refine_Polyhedron_3.h>
#include <CGAL/internal/Hole_filling/Triangulate_hole_Polyhedron_3.h>
#include <CGAL/Default.h>
#include <vector>
#include <boost/tuple/tuple.hpp>

namespace CGAL {

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*!
\ingroup PkgPolygonMeshProcessing
@brief Function triangulating and refining a hole in a polygon mesh.

@tparam PolygonMesh must be model of `MutableFaceGraph`
@tparam FacetOutputIterator iterator holding `boost::graph_traits<PolygonMesh>::face_descriptor` for patch faces.
@tparam VertexOutputIterator iterator holding `boost::graph_traits<PolygonMesh>::vertex_descriptor` for patch vertices.

@param pmesh polygon mesh which has the hole
@param border_halfedge a border halfedge incident to the hole
@param face_out iterator over patch facess
@param vertex_out iterator over patch vertices without including the boundary
@param density_control_factor factor for density where larger values cause denser refinements
@param use_delaunay_triangulation if `true`, use the Delaunay triangulation face search space

@return pair of @a face_out and @a vertex_out

\todo handle islands
\todo use Delaunay by default
*/
template<
  class PolygonMesh,
  class FaceOutputIterator,
  class VertexOutputIterator
>
std::pair<FaceOutputIterator, VertexOutputIterator> 
triangulate_and_refine_hole(PolygonMesh& pmesh, 
                            typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge, 
                            FaceOutputIterator face_out,
                            VertexOutputIterator vertex_out,
                            double density_control_factor = std::sqrt(2.0),
                            bool use_delaunay_triangulation = false) 
{
  std::vector<typename boost::graph_traits<PolygonMesh>::face_descriptor> patch;
  triangulate_hole(pmesh, border_halfedge, std::back_inserter(patch), use_delaunay_triangulation);
  face_out = std::copy(patch.begin(), patch.end(), face_out);
  return refine(pmesh, patch.begin(), patch.end(), face_out, vertex_out, density_control_factor);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
\ingroup PkgPolygonMeshProcessing
@brief Function triangulating, refining and fairing a hole in a polygon mesh.

If \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available
and `CGAL_EIGEN3_ENABLED` is defined, an overload of this function is available
with `SparseLinearSolver` being:
\code
  CGAL::Eigen_solver_traits<
  Eigen::SparseLU<
  CGAL::Eigen_sparse_matrix<double>::EigenType,
  Eigen::COLAMDOrdering<int> >  >
\endcode
and `WeightCalculator` being `CGAL::internal::Cotangent_weight_with_voronoi_area_fairing<PolygonMesh>`.
For using an alternative model of `FairWeightCalculator` with the default solver,
one can pass `CGAL::Default()` as `solver`.

@tparam SparseLinearSolver a model of `SparseLinearAlgebraTraitsWithPreFactor_d`
@tparam WeightCalculator a model of `FairWeightCalculator`
@tparam PolygonMesh must be model of  `MutableFaceGraph`
@tparam FaceOutputIterator iterator holding `boost::graph_traits<PolygonMesh>::face_descriptor` for patch faces.
@tparam VertexOutputIterator iterator holding `boost::graph_traits<PolygonMesh>::vertex_descriptor` for patch vertices.

@param pmesh a polygon mesh which has the hole
@param border_halfedge a border halfedge incident to the hole
@param face_out iterator over patch faces
@param vertex_out iterator over patch vertices without including the boundary
@param weight_calculator function object to calculate weights, default to Cotangent weights and can be omitted
@param density_control_factor factor for density where larger values cause denser refinements
@param continuity tangential continuity, default to `FAIRING_C_1` and can be omitted
@param use_delaunay_triangulation if `true`, use the Delaunay triangulation face search space
@param solver An instance of the sparse linear solver to use. Note that the current implementation is
              not using the value passed but the default constructed one.

@return tuple of 
 - bool: `true` if fairing is successful
 - @a face_out
 - @a vertex_out
 
\todo handle islands

\todo use Delaunay by default
\todo WeightCalculator should be a property map
 */
template<
  class WeightCalculator,
  class SparseLinearSolver,
  class PolygonMesh,
  class FaceOutputIterator,
  class VertexOutputIterator
>
boost::tuple<bool, FaceOutputIterator, VertexOutputIterator>
triangulate_refine_and_fair_hole(PolygonMesh& pmesh, 
                                 typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge, 
                                 FaceOutputIterator face_out,
                                 VertexOutputIterator vertex_out,
                                 WeightCalculator weight_calculator,
                                 SparseLinearSolver
                                 #ifdef DOXYGEN_RUNNING
                                 solver,
                                 #else
                                 /* solver */,
                                 #endif
                                 double density_control_factor = std::sqrt(2.0),
                                 bool use_delaunay_triangulation = false,
                                 Fairing_continuity continuity = FAIRING_C_1)
{
  std::vector<typename boost::graph_traits<PolygonMesh>::vertex_descriptor> patch;

  face_out = triangulate_and_refine_hole
    (pmesh, border_halfedge, face_out, std::back_inserter(patch), density_control_factor, use_delaunay_triangulation)
              .first;

  typedef CGAL::internal::Fair_default_sparse_linear_solver::Solver Default_solver;
  typedef typename Default::Get<SparseLinearSolver, Default_solver>::type Solver;

  bool fair_success = fair<Solver>(pmesh, patch.begin(), patch.end(), weight_calculator, continuity);

  vertex_out = std::copy(patch.begin(), patch.end(), vertex_out);
  return boost::make_tuple(fair_success, face_out, vertex_out);
}

//use default SparseLinearSolver and WeightCalculator
template<
  class PolygonMesh,
  class FaceOutputIterator,
  class VertexOutputIterator
>
boost::tuple<bool, FaceOutputIterator, VertexOutputIterator>
triangulate_refine_and_fair_hole(PolygonMesh& pmesh, 
                                 typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge, 
                                 FaceOutputIterator face_out,
                                 VertexOutputIterator vertex_out,
                                 double density_control_factor = std::sqrt(2.0),
                                 bool use_delaunay_triangulation = false,
                                 Fairing_continuity continuity = FAIRING_C_1)
{
  CGAL::internal::Cotangent_weight_with_voronoi_area_fairing<PolygonMesh> wc(pmesh);

  return triangulate_refine_and_fair_hole
    (pmesh, border_halfedge, face_out, vertex_out, wc, Default(), density_control_factor, use_delaunay_triangulation, continuity);
}


} // namespace CGAL

#endif
