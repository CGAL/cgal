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
@brief Function triangulating and refining a hole in a surface mesh.

@tparam Polyhedron a \cgal polyhedron
@tparam FacetOutputIterator iterator holding `Polyhedron::Facet_handle` for patch facets.
@tparam VertexOutputIterator iterator holding `Polyhedron::Vertex_handle` for patch vertices.

@param polyhedron surface mesh which has the hole
@param border_halfedge a border halfedge incident to the hole
@param facet_out iterator over patch facets
@param vertex_out iterator over patch vertices without including the boundary
@param density_control_factor factor for density where larger values cause denser refinements
@param use_delaunay_triangulation if `true`, use the Delaunay triangulation facet search space

@return pair of @a facet_out and @a vertex_out

\todo handle islands
\todo `Polyhedron` should be a model of `MutableFaceGraph`
\todo use Delaunay by default
*/
template<
  class Polyhedron,
  class FacetOutputIterator,
  class VertexOutputIterator
>
std::pair<FacetOutputIterator, VertexOutputIterator> 
triangulate_and_refine_hole(Polyhedron& polyhedron, 
                            typename boost::graph_traits<Polyhedron>::halfedge_descriptor border_halfedge, 
                            FacetOutputIterator facet_out,
                            VertexOutputIterator vertex_out,
                            double density_control_factor = std::sqrt(2.0),
                            bool use_delaunay_triangulation = false) 
{
  std::vector<typename boost::graph_traits<Polyhedron>::face_descriptor> patch;
  triangulate_hole(polyhedron, border_halfedge, std::back_inserter(patch), use_delaunay_triangulation);
  facet_out = std::copy(patch.begin(), patch.end(), facet_out);
  return refine(polyhedron, patch.begin(), patch.end(), facet_out, vertex_out, density_control_factor);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
\ingroup PkgPolygonMeshProcessing
@brief Function triangulating, refining and fairing a hole in surface mesh.

If \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available
and `CGAL_EIGEN3_ENABLED` is defined, an overload of this function is available
with `SparseLinearSolver` being:
\code
  CGAL::Eigen_solver_traits<
  Eigen::SparseLU<
  CGAL::Eigen_sparse_matrix<double>::EigenType,
  Eigen::COLAMDOrdering<int> >  >
\endcode
and `WeightCalculator` being `CGAL::internal::Cotangent_weight_with_voronoi_area_fairing<Polyhedron>`.
For using an alternative model of `FairWeightCalculator` with the default solver,
one can pass `CGAL::Default()` as `solver`.

@tparam SparseLinearSolver a model of `SparseLinearAlgebraTraitsWithPreFactor_d`
@tparam WeightCalculator a model of `FairWeightCalculator`
@tparam Polyhedron a \cgal polyhedron
@tparam FacetOutputIterator iterator holding `Polyhedron::Facet_handle` for patch facets.
@tparam VertexOutputIterator iterator holding `Polyhedron::Vertex_handle` for patch vertices.

@param polyhedron surface mesh which has the hole
@param border_halfedge a border halfedge incident to the hole
@param facet_out iterator over patch facets
@param vertex_out iterator over patch vertices without including the boundary
@param weight_calculator function object to calculate weights, default to Cotangent weights and can be omitted
@param density_control_factor factor for density where larger values cause denser refinements
@param continuity tangential continuity, default to `FAIRING_C_1` and can be omitted
@param use_delaunay_triangulation if `true`, use the Delaunay triangulation facet search space
@param solver An instance of the sparse linear solver to use. Note that the current implementation is
              not using the value passed but the default constructed one.

@return tuple of 
 - bool: `true` if fairing is successful
 - @a facet_out
 - @a vertex_out
 
\todo handle islands
\todo `Polyhedron` should be a model of `MutableFaceGraph`
\todo use Delaunay by default
\todo WeightCalculator should be a property map
 */
template<
  class WeightCalculator,
  class SparseLinearSolver,
  class Polyhedron,
  class FacetOutputIterator,
  class VertexOutputIterator
>
boost::tuple<bool, FacetOutputIterator, VertexOutputIterator>
triangulate_refine_and_fair_hole(Polyhedron& polyhedron, 
                                 typename boost::graph_traits<Polyhedron>::halfedge_descriptor border_halfedge, 
                                 FacetOutputIterator facet_out,
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
  std::vector<typename boost::graph_traits<Polyhedron>::vertex_descriptor> patch;

  facet_out = triangulate_and_refine_hole
    (polyhedron, border_halfedge, facet_out, std::back_inserter(patch), density_control_factor, use_delaunay_triangulation)
              .first;

  typedef CGAL::internal::Fair_default_sparse_linear_solver::Solver Default_solver;
  typedef typename Default::Get<SparseLinearSolver, Default_solver>::type Solver;

  bool fair_success = fair<Solver>(polyhedron, patch.begin(), patch.end(), weight_calculator, continuity);

  vertex_out = std::copy(patch.begin(), patch.end(), vertex_out);
  return boost::make_tuple(fair_success, facet_out, vertex_out);
}

//use default SparseLinearSolver and WeightCalculator
template<
  class Polyhedron,
  class FacetOutputIterator,
  class VertexOutputIterator
>
boost::tuple<bool, FacetOutputIterator, VertexOutputIterator>
triangulate_refine_and_fair_hole(Polyhedron& polyhedron, 
                                 typename boost::graph_traits<Polyhedron>::halfedge_descriptor border_halfedge, 
                                 FacetOutputIterator facet_out,
                                 VertexOutputIterator vertex_out,
                                 double density_control_factor = std::sqrt(2.0),
                                 bool use_delaunay_triangulation = false,
                                 Fairing_continuity continuity = FAIRING_C_1)
{
  CGAL::internal::Cotangent_weight_with_voronoi_area_fairing<Polyhedron> wc;

  return triangulate_refine_and_fair_hole
    (polyhedron, border_halfedge, facet_out, vertex_out, wc, Default(), density_control_factor, use_delaunay_triangulation, continuity);
}


} // namespace CGAL

#endif
