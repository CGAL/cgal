#ifndef CGAL_FILL_HOLE_POLYHEDRON_H
#define CGAL_FILL_HOLE_POLYHEDRON_H

// Helper functions which combine triangulate, fair, and refine 

#include <CGAL/internal/Fair_Polyhedron_3.h>
#include <CGAL/internal/Refine_Polyhedron_3.h>
#include <CGAL/internal/Triangulate_hole_Polyhedron_3.h>
#include <vector>

namespace CGAL {

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Function triangulating and refining a hole in surface mesh.
 * 
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam FacetOutputIterator iterator holding 'Polyhedron::Facet_handle' for patch facets.
 * @tparam VertexOutputIterator iterator holding 'Polyhedron::Vertex_handle' for patch vertices.
 *
 * @param[in, out] polyhedron surface mesh which has the hole
 * @param border_halfedge a border halfedge incident to the hole
 * @param[out] facet_out iterator over patch facets
 * @param[out] vertex_out iterator over patch vertices without including boundary
 * @param density_control_factor factor for density where larger values cause denser refinements
 * 
 * \warning Using this function on very large holes might not be feasible, since the cost of triangulation is O(n^3).
 */
template<class Polyhedron, class FacetOutputIterator, class VertexOutputIterator>
void triangulate_and_refine_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  FacetOutputIterator facet_out,
  VertexOutputIterator vertex_out,
  double density_control_factor = std::sqrt(2.0)) 
{
  std::vector<typename Polyhedron::Facet_handle> patch;
  triangulate_hole(polyhedron, border_halfedge, std::back_inserter(patch));

  refine(polyhedron, patch.begin(), patch.end(), facet_out, vertex_out, density_control_factor);
}

//// do not use OutputIterator
//template<class Polyhedron>
//void triangulate_and_refine_hole(Polyhedron& polyhedron, 
//  typename Polyhedron::Halfedge_handle border_halfedge, 
//  double density_control_factor = std::sqrt(2.0)
//  ) 
//{
//  triangulate_and_refine_hole(polyhedron, border_halfedge,
//    boost::make_function_output_iterator(CGAL::internal::Nop_functor()),
//    density_control_factor);
//}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** 
 * @brief Function triangulating, refining and fairing a hole in surface mesh.
 *
 * @tparam SparseLinearSolver a model of SparseLinearAlgebraTraitsWithPreFactor_d and can be omitted if Eigen defined...(give exact models etc)
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam WeightCalculator a model of "weight model" and default to Cotangent weights
 *
 * @param polyhedron surface mesh which has the hole
 * @param border_halfedge a border halfedge incident to the hole
 * @param[out] facet_out iterator over patch facets
 * @param[out] vertex_out iterator over patch vertices without including boundary
 * @param density_control_factor factor for density where larger values cause denser refinements
 * @param weight_calculator function object to calculate weights
 */
template<class SparseLinearSolver, class WeightCalculator, class Polyhedron, class FacetOutputIterator, class VertexOutputIterator>
void triangulate_refine_and_fair_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  FacetOutputIterator facet_out,
  VertexOutputIterator vertex_out,
  WeightCalculator weight_calculator,
  double density_control_factor = std::sqrt(2.0)
  )
{
  std::vector<typename Polyhedron::Vertex_handle> patch;
  triangulate_and_refine_hole(polyhedron, border_halfedge, facet_out, std::back_inserter(patch), density_control_factor);
  fair<SparseLinearSolver>(polyhedron, patch.begin(), patch.end(), weight_calculator);
  std::copy(patch.begin(), patch.end(), vertex_out);
}

//use default SparseLinearSolver
template<class WeightCalculator, class Polyhedron, class FacetOutputIterator, class VertexOutputIterator>
void triangulate_refine_and_fair_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  FacetOutputIterator facet_out,
  VertexOutputIterator vertex_out,
  WeightCalculator weight_calculator,
  double density_control_factor = std::sqrt(2.0)
  )
{
  typedef CGAL::internal::Fair_default_sparse_linear_solver::Solver Sparse_linear_solver;
  triangulate_refine_and_fair_hole<Sparse_linear_solver, WeightCalculator, Polyhedron, FacetOutputIterator, VertexOutputIterator>
    (polyhedron, border_halfedge, facet_out, vertex_out, weight_calculator, density_control_factor);
}

//use default WeightCalculator
template<class SparseLinearSolver, class Polyhedron, class FacetOutputIterator, class VertexOutputIterator>
void triangulate_refine_and_fair_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  FacetOutputIterator facet_out,
  VertexOutputIterator vertex_out,
  double density_control_factor = std::sqrt(2.0)
  )
{
  typedef CGAL::Fairing_cotangent_weight<Polyhedron> Weight_calculator;
  triangulate_refine_and_fair_hole<SparseLinearSolver, Weight_calculator, Polyhedron, FacetOutputIterator, VertexOutputIterator>
    (polyhedron, border_halfedge, facet_out, vertex_out,  Weight_calculator(), density_control_factor);
}

//use default SparseLinearSolver and WeightCalculator
template<class Polyhedron, class FacetOutputIterator, class VertexOutputIterator>
void triangulate_refine_and_fair_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  FacetOutputIterator facet_out,
  VertexOutputIterator vertex_out,
  double density_control_factor = std::sqrt(2.0)
  )
{
  typedef CGAL::internal::Fair_default_sparse_linear_solver::Solver Sparse_linear_solver;
  triangulate_refine_and_fair_hole<Sparse_linear_solver, Polyhedron, FacetOutputIterator, VertexOutputIterator>
    (polyhedron, border_halfedge, facet_out, vertex_out, density_control_factor);
}

} // namespace CGAL

#endif