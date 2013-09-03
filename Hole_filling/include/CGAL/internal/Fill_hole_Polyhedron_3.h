#ifndef CGAL_HOLE_FILLING_FILL_HOLE_POLYHEDRON_H
#define CGAL_HOLE_FILLING_FILL_HOLE_POLYHEDRON_H

// Helper functions which combine triangulate, fair, and refine 

#include <CGAL/internal/Fair_Polyhedron_3.h>
#include <CGAL/internal/Refine_Polyhedron_3.h>
#include <CGAL/internal/Triangulate_hole_Polyhedron_3.h>
#include <vector>
#include <boost/tuple/tuple.hpp>

namespace CGAL {

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*!
\ingroup PkgHoleFilling
@brief Function triangulating and refining a hole in surface mesh.

@tparam Polyhedron a %CGAL polyhedron
@tparam FacetOutputIterator iterator holding `Polyhedron::Facet_handle` for patch facets.
@tparam VertexOutputIterator iterator holding `Polyhedron::Vertex_handle` for patch vertices.

@param polyhedron surface mesh which has the hole
@param border_halfedge a border halfedge incident to the hole
@param facet_out iterator over patch facets
@param vertex_out iterator over patch vertices without including boundary
@param density_control_factor factor for density where larger values cause denser refinements

@return pair of @a facet_out and @a vertex_out
*/
template<class Polyhedron, class FacetOutputIterator, class VertexOutputIterator>
std::pair<FacetOutputIterator, VertexOutputIterator> 
triangulate_and_refine_hole(
  Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  FacetOutputIterator facet_out,
  VertexOutputIterator vertex_out,
  double density_control_factor = std::sqrt(2.0),
  bool use_delaunay_triangulation = false
  ) 
{
  std::vector<typename Polyhedron::Facet_handle> patch;
  triangulate_hole(polyhedron, border_halfedge, std::back_inserter(patch), use_delaunay_triangulation);
  facet_out = std::copy(patch.begin(), patch.end(), facet_out);
  return refine(polyhedron, patch.begin(), patch.end(), facet_out, vertex_out, density_control_factor);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
\ingroup PkgHoleFilling
@brief Function triangulating, refining and fairing a hole in surface mesh.

@tparam SparseLinearSolver a model of SparseLinearAlgebraTraitsWithPreFactor_d and can be omitted if Eigen defined...(give exact models etc)
@tparam WeightCalculator a model of FairWeightCalculator and can be omitted to use default Cotangent weights
@tparam Polyhedron a %CGAL polyhedron
@tparam FacetOutputIterator iterator holding `Polyhedron::Facet_handle` for patch facets.
@tparam VertexOutputIterator iterator holding `Polyhedron::Vertex_handle` for patch vertices.

@param polyhedron surface mesh which has the hole
@param border_halfedge a border halfedge incident to the hole
@param facet_out iterator over patch facets
@param vertex_out iterator over patch vertices without including boundary
@param weight_calculator function object to calculate weights, default to Cotangent weights and can be omitted
@param density_control_factor factor for density where larger values cause denser refinements
@param continuity tangential continuity, default to FAIRING_C_1 and can be omitted
@return tuple of 
 - bool: true if fairing is successful
 - @a facet_out
 - @a vertex_out
 */
template<class SparseLinearSolver, class WeightCalculator, class Polyhedron, class FacetOutputIterator, class VertexOutputIterator>
boost::tuple<bool, FacetOutputIterator, VertexOutputIterator>
triangulate_refine_and_fair_hole(
  Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  FacetOutputIterator facet_out,
  VertexOutputIterator vertex_out,
  WeightCalculator weight_calculator,
  double density_control_factor = std::sqrt(2.0),
  bool use_delaunay_triangulation = false,
  Fairing_continuity continuity = FAIRING_C_1
  )
{
  std::vector<typename Polyhedron::Vertex_handle> patch;

  facet_out = triangulate_and_refine_hole
    (polyhedron, border_halfedge, facet_out, std::back_inserter(patch), density_control_factor, use_delaunay_triangulation)
              .first;

  bool fair_success = fair<SparseLinearSolver>(polyhedron, patch.begin(), patch.end(), weight_calculator, continuity);

  vertex_out = std::copy(patch.begin(), patch.end(), vertex_out);
  return boost::make_tuple(fair_success, facet_out, vertex_out);
}

//use default SparseLinearSolver
template<class WeightCalculator, class Polyhedron, class FacetOutputIterator, class VertexOutputIterator>
boost::tuple<bool, FacetOutputIterator, VertexOutputIterator>
triangulate_refine_and_fair_hole(
  Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  FacetOutputIterator facet_out,
  VertexOutputIterator vertex_out,
  WeightCalculator weight_calculator,
  double density_control_factor = std::sqrt(2.0),
  bool use_delaunay_triangulation = false,
  Fairing_continuity continuity = FAIRING_C_1
  )
{
  typedef CGAL::internal::Fair_default_sparse_linear_solver::Solver Sparse_linear_solver;
  return triangulate_refine_and_fair_hole<Sparse_linear_solver, WeightCalculator, Polyhedron, FacetOutputIterator, VertexOutputIterator>
    (polyhedron, border_halfedge, facet_out, vertex_out, weight_calculator, density_control_factor, use_delaunay_triangulation, continuity);
}

//use default WeightCalculator
template<class SparseLinearSolver, class Polyhedron, class FacetOutputIterator, class VertexOutputIterator>
boost::tuple<bool, FacetOutputIterator, VertexOutputIterator>
triangulate_refine_and_fair_hole(
  Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  FacetOutputIterator facet_out,
  VertexOutputIterator vertex_out,
  double density_control_factor = std::sqrt(2.0),
  bool use_delaunay_triangulation = false,
  Fairing_continuity continuity = FAIRING_C_1
  )
{
  typedef CGAL::internal::Cotangent_weight_with_voronoi_area_fairing<Polyhedron> Weight_calculator;
  return triangulate_refine_and_fair_hole<SparseLinearSolver, Weight_calculator, Polyhedron, FacetOutputIterator, VertexOutputIterator>
    (polyhedron, border_halfedge, facet_out, vertex_out,  Weight_calculator(), density_control_factor, use_delaunay_triangulation, continuity);
}

//use default SparseLinearSolver and WeightCalculator
template<class Polyhedron, class FacetOutputIterator, class VertexOutputIterator>
boost::tuple<bool, FacetOutputIterator, VertexOutputIterator>
triangulate_refine_and_fair_hole(
  Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  FacetOutputIterator facet_out,
  VertexOutputIterator vertex_out,
  double density_control_factor = std::sqrt(2.0),
  bool use_delaunay_triangulation = false,
  Fairing_continuity continuity = FAIRING_C_1
  )
{
  typedef CGAL::internal::Fair_default_sparse_linear_solver::Solver Sparse_linear_solver;
  return triangulate_refine_and_fair_hole<Sparse_linear_solver, Polyhedron, FacetOutputIterator, VertexOutputIterator>
    (polyhedron, border_halfedge, facet_out, vertex_out, density_control_factor, use_delaunay_triangulation, continuity);
}

} // namespace CGAL

#endif