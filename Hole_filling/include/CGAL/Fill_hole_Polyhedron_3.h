// adapted from fill.cpp -IOY
#ifndef CGAL_FILL_HOLE_POLYHEDRON_H
#define CGAL_FILL_HOLE_POLYHEDRON_H
#include <CGAL/internal/Fair.h>
#include <CGAL/internal/Refine.h>
#include <CGAL/internal/Triangulate.h>

#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <CGAL/internal/Weights.h>
#include <CGAL/trace.h>

#include <vector>
#include <limits>
#include <set>
#include <map>

namespace CGAL {
namespace internal {

  template<class Polyhedron, class SparseLinearSolver, class WeightCalculator>
  class Fill_hole_Polyhedron_3 {

    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Facet_handle Facet_handle;
    typedef typename Polyhedron::Halfedge_around_facet_circulator  Halfedge_around_facet_circulator;
    typedef typename Polyhedron::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;

    //members
    bool is_refining;
    bool is_fairing;
    double alpha;
    WeightCalculator weight_calculator;

  public:
    Fill_hole_Polyhedron_3(bool is_refining, double alpha, bool is_fairing, WeightCalculator weight_calculator = WeightCalculator()) 
      : is_refining(is_refining), is_fairing(is_fairing), alpha(alpha), weight_calculator(weight_calculator)
    { }

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

  public:
    void operator()(Polyhedron& poly, Halfedge_handle it)
    {
      if(! it->is_border()){
        return;
      }

      // save boundary vertices and scale attributes to use in refining and fairing state
      std::map<Vertex_handle, double> scale_attribute;
      std::set<Vertex_handle> boundary_vertices;
      if(is_refining) {
        Halfedge_around_facet_circulator circ(it), done(circ);
        do{
          average_length(poly, circ->vertex(), scale_attribute); 
          if(is_fairing){ boundary_vertices.insert(circ->vertex()); }
        } while (++circ != done);
      }

      // Triangulate
      std::set<Facet_handle> facets;
      internal::Triangulate_Hole_Polyhedron_3<Polyhedron> triangulate_functor;
      triangulate_functor(poly, it, facets);

      if(is_refining) {
        // Refine
        internal::Refine_Polyhedron_3<Polyhedron> refine_functor(alpha);
        refine_functor(poly, scale_attribute, facets);
        CGAL_TRACE_STREAM << "after refine |facets| = " << facets.size() << std::endl;

        if(is_fairing) {
          std::set<Vertex_handle> interior_vertices;
          for(std::set<Facet_handle>::iterator it = facets.begin(); it != facets.end(); ++it) {
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
          internal::Fair_Polyhedron_3<Polyhedron, SparseLinearSolver, WeightCalculator> fair_functor(weight_calculator);
          fair_functor(poly, interior_vertices);
        } // if(is_fairing)
      } // if(is_refining)
    }
  };
}//namespace internal

  template<class SparseLinearSolver, class Polyhedron, class WeightCalculator>
  void fill(Polyhedron& poly, 
    typename Polyhedron::Halfedge_handle it, 
    bool refine = true,
    double density_control_factor = 1.41 /* ~sqrt(2) */,
    bool fair = true,
    WeightCalculator weight_calculator = WeightCalculator()
    )
  {
    internal::Fill_hole_Polyhedron_3<Polyhedron, SparseLinearSolver, WeightCalculator> 
      fill_functor(refine, density_control_factor, fair, weight_calculator);
    fill_functor(poly, it);
  }

  template<class SparseLinearSolver, class Polyhedron>
  void fill(Polyhedron& poly, 
    typename Polyhedron::Halfedge_handle it, 
    bool refine = true,
    double density_control_factor = 1.41 /* ~sqrt(2) */,
    bool fair = true
    )
  {
    typedef typename CGAL::internal::Fairing_weight_selector<Polyhedron, COTANGENT_WEIGHTING>::weight_calculator
      Weight_calculator;
    fill<SparseLinearSolver, Polyhedron, Weight_calculator>
      (poly, it, refine, density_control_factor, fair, Weight_calculator());
  }

  template<class Polyhedron, class WeightCalculator>
  void fill(Polyhedron& poly, 
    typename Polyhedron::Halfedge_handle it, 
    bool refine = true,
    double density_control_factor = 1.41 /* ~sqrt(2) */,
    bool fair = true,
    WeightCalculator weight_calculator = WeightCalculator()
    )
  {
    typedef CGAL::internal::Fair_default_sparse_linear_solver::Solver Sparse_linear_solver;
    fill<Sparse_linear_solver, Polyhedron, WeightCalculator>
      (poly, it, refine, density_control_factor, fair, weight_calculator);
  }

  template<class Polyhedron>
  void fill(Polyhedron& poly, 
    typename Polyhedron::Halfedge_handle it, 
    bool refine = true,
    double density_control_factor = 1.41 /* ~sqrt(2) */,
    bool fair = true
    )
  {
    typedef CGAL::internal::Fair_default_sparse_linear_solver::Solver Sparse_linear_solver;
    fill<Sparse_linear_solver, Polyhedron>
      (poly, it, refine, density_control_factor, fair);
  }
} // namespace CGAL

#endif