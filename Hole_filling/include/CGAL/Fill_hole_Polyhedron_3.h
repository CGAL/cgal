// adapted from fill.cpp -IOY
#ifndef CGAL_FILL_HOLE_POLYHEDRON_H
#define CGAL_FILL_HOLE_POLYHEDRON_H
#include <CGAL/internal/Fair.h>
#include <CGAL/internal/Refine.h>
#include <CGAL/internal/Triangulate.h>

#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <CGAL/internal/Weights.h>
#include <CGAL/trace.h>

#include <boost/function_output_iterator.hpp>

#include <vector>
#include <limits>
#include <set>
#include <map>

namespace CGAL {
namespace internal {

  struct Nop_functor {
    template<class T>
    void operator()(const T& t) const {}
  };

  template<class Polyhedron>
  class Fill_hole_Polyhedron_3 {

    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Facet_handle Facet_handle;
    typedef typename Polyhedron::Halfedge_around_facet_circulator  Halfedge_around_facet_circulator;
    typedef typename Polyhedron::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;

  private:
    void average_length(Polyhedron& poly, Vertex_handle vh, std::map<Vertex_handle, double>& scale_attribute)
    {
      const Point_3& vp = vh->point(); 
      Halfedge_around_vertex_circulator circ(vh->vertex_begin()), done(circ);
      int deg = 0;
      double sum = 0;
      do {
        const Point_3& vq = circ->opposite()->vertex()->point();
        sum += std::sqrt(CGAL::squared_distance(vp, vq));
        ++deg;
        ++circ;
      } while(circ != done);
      scale_attribute[vh] = sum/deg;
    }

    void get_boundary_vertices(Polyhedron& poly, Halfedge_handle it, std::set<Vertex_handle>& boundary_vertices) {
      Halfedge_around_facet_circulator circ(it), done(circ);
      do{
        boundary_vertices.insert(circ->vertex()); 
      } while (++circ != done);
    }

  public:
    template<class OutputIterator>
    void triangulate_and_refine_hole(Polyhedron& poly, Halfedge_handle it, double alpha, OutputIterator output)
    {
      if(! it->is_border()){ return; }

      // compute scale_attribute before triangulation
      std::map<Vertex_handle, double> scale_attribute;
      Halfedge_around_facet_circulator circ(it), done(circ);
      do{
        average_length(poly, circ->vertex(), scale_attribute); 
      } while (++circ != done);
      // Triangulate
      std::set<Facet_handle> facets;
      triangulate_hole(poly, it, std::inserter(facets, facets.begin()));
      // Refine
      internal::Refine_Polyhedron_3<Polyhedron> refine_functor(alpha);
      refine_functor(poly, scale_attribute, facets);

      std::copy(facets.begin(), facets.end(), output);
    }

    template<class SparseLinearSolver, class WeightCalculator>
    void triangulate_refine_and_fair
      (Polyhedron& poly, Halfedge_handle it, double alpha, WeightCalculator weight_calculator) 
    {
      if(! it->is_border()){ return; }

      // save boundary vertices before triangulation
      std::set<Vertex_handle> boundary_vertices;
      get_boundary_vertices(poly, it, boundary_vertices);

      //Triangulate and refine
      std::vector<Facet_handle> facets;
      triangulate_and_refine_hole(poly, it, alpha, std::back_inserter(facets));

      // get interior vertices
      std::set<Vertex_handle> interior_vertices;
      for(std::vector<Facet_handle>::iterator it = facets.begin(); it != facets.end(); ++it) {
        Halfedge_around_facet_circulator circ = (*it)->facet_begin();
        do {
          if(boundary_vertices.find(circ->vertex()) == boundary_vertices.end()) {
            interior_vertices.insert(circ->vertex());
          }
        } while(++circ != (*it)->facet_begin());
      }

      CGAL_TRACE_STREAM << "before fair |boundary vertices| = " << boundary_vertices.size() << std::endl;
      CGAL_TRACE_STREAM << "before fair |interior vertices| = " << interior_vertices.size() << std::endl;
      // Fair
      fair<SparseLinearSolver>(poly, interior_vertices, weight_calculator);
    }
  };
}//namespace internal

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
  internal::Triangulate_Hole_Polyhedron_3<Polyhedron>()(polyhedron, border_halfedge, output);
}

template<class Polyhedron>
void triangulate_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge
  ) 
{
  internal::Triangulate_Hole_Polyhedron_3<Polyhedron>()(polyhedron, border_halfedge);
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
  fill_functor.triangulate_refine_and_fair<Sparse_linear_solver, WeightCalculator>
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
  fill_functor.triangulate_refine_and_fair<SparseLinearSolver, WeightCalculator>
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