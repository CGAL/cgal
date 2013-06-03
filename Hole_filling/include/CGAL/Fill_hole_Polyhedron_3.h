#ifndef CGAL_FILL_HOLE_POLYHEDRON_H
#define CGAL_FILL_HOLE_POLYHEDRON_H
// This file is the place where all free functions related to hole triangulation stand
#include <CGAL/internal/Fill_hole_connector.h>

namespace CGAL {

/**
 * @brief Function triangulating a hole in surface mesh.
 * 
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam OutputIterator iterator holding 'Polyhedron::Facet_handle' for patch facets.
 * @param[in, out] polyhedron surface mesh which has the hole
 * @param border_halfedge a border halfedge incident to the hole
 * @param[out] output iterator over patch facets. It can be omitted if output is not required.
 * 
 * \warning Using this function on very large holes might not be feasible, since the cost of triangulation is O(n^3).
 */
template<class Polyhedron, class OutputIterator>
void triangulate_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  OutputIterator output
  ) 
{
  internal::Triangulate_hole_Polyhedron_3<Polyhedron>()(polyhedron, border_halfedge, output);
}

template<class Polyhedron>
void triangulate_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge
  ) 
{
  triangulate_hole(polyhedron, border_halfedge, 
    boost::make_function_output_iterator(CGAL::internal::Nop_functor()));
}

/**
 * @brief Function triangulating and refining a hole in surface mesh.
 * 
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam OutputIterator iterator holding 'Polyhedron::Facet_handle' for patch facets.
 * @param[in, out] polyhedron surface mesh which has the hole
 * @param border_halfedge a border halfedge incident to the hole
 * @param[out] output iterator over patch facets. It can be omitted if output is not required.
 * @param density_control_factor factor for density where larger values cause denser refinements
 * 
 * \warning Using this function on very large holes might not be feasible, since the cost of triangulation is O(n^3).
 */
template<class Polyhedron, class OutputIterator>
void triangulate_and_refine_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  OutputIterator output,
  double density_control_factor = std::sqrt(2.0)) 
{
  internal::Fill_hole_Polyhedron_3<Polyhedron>().
    triangulate_and_refine_hole(polyhedron, border_halfedge, density_control_factor, output);
}

template<class Polyhedron>
void triangulate_and_refine_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  double density_control_factor = std::sqrt(2.0)
  ) 
{
  triangulate_and_refine_hole(polyhedron, border_halfedge,
    boost::make_function_output_iterator(CGAL::internal::Nop_functor()),
    density_control_factor);
}

//////////////////////////////////////////////////////////////////////////

template<class Polyhedron, class WeightCalculator>
void triangulate_refine_and_fair_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  double density_control_factor = std::sqrt(2.0),
  WeightCalculator weight_calculator = WeightCalculator()
  )
{
  typedef CGAL::internal::Fair_default_sparse_linear_solver::Solver Sparse_linear_solver;
  //triangulate_refine_and_fair_hole<Sparse_linear_solver, Polyhedron, WeightCalculator>
  //  (polyhedron, border_halfedge, new_facets, density_control_factor, weight_calculator);

  CGAL::internal::Fill_hole_Polyhedron_3<Polyhedron> fill_functor;
  fill_functor.template triangulate_refine_and_fair<Sparse_linear_solver, WeightCalculator>
    (polyhedron, border_halfedge, density_control_factor, weight_calculator);
}

/** 
 * @brief Function triangulating, refining and fairing a hole in surface mesh.
 *
 * @tparam SparseLinearSolver a model of SparseLinearAlgebraTraitsWithPreFactor_d and can be omitted if Eigen defined...(give exact models etc)
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam WeightCalculator a model of "weight model" and default to Cotangent weights
 * @param polyhedron surface mesh which has the hole
 * @param border_halfedge a border halfedge incident to the hole
 * @param[out] new_facets patch facets to close the hole
 * @param density_control_factor factor for density where larger values cause denser refinements
 */
template<class SparseLinearSolver, class Polyhedron, class WeightCalculator>
void triangulate_refine_and_fair_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  double density_control_factor = std::sqrt(2.0),
  WeightCalculator weight_calculator = WeightCalculator()
  )
{
  CGAL::internal::Fill_hole_Polyhedron_3<Polyhedron> fill_functor;
  fill_functor.template triangulate_refine_and_fair<SparseLinearSolver, WeightCalculator>
    (polyhedron, border_halfedge, new_facets, density_control_factor, weight_calculator);
}

//template<class SparseLinearSolver, class Polyhedron>
//void triangulate_refine_and_fair_hole(Polyhedron& polyhedron, 
//  typename Polyhedron::Halfedge_handle border_halfedge, 
//  std::set<typename Polyhedron::Facet_handle>* new_facets = NULL,
//  double density_control_factor = std::sqrt(2.0)
//  )
//{
//  typedef typename CGAL::internal::Fairing_weight_selector<Polyhedron, COTANGENT_WEIGHTING>::weight_calculator
//    Weight_calculator;
//  triangulate_refine_and_fair_hole<SparseLinearSolver, Polyhedron, Weight_calculator>
//    (polyhedron, border_halfedge, new_facets, density_control_factor, Weight_calculator());
//}

//template<class Polyhedron>
//void triangulate_refine_and_fair_hole(Polyhedron& polyhedron, 
//  typename Polyhedron::Halfedge_handle border_halfedge, 
//  std::set<typename Polyhedron::Facet_handle>* new_facets = NULL,
//  double density_control_factor = std::sqrt(2.0)
//  )
//{
//  typedef CGAL::internal::Fair_default_sparse_linear_solver::Solver Sparse_linear_solver;
//  triangulate_refine_and_fair_hole<Sparse_linear_solver, Polyhedron>
//    (polyhedron, border_halfedge, new_facets, density_control_factor);
//}
} // namespace CGAL

#endif