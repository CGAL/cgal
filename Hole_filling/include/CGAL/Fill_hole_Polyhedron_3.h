// adapted from fill.cpp -IOY
#ifndef CGAL_FILL_HOLE_POLYHEDRON_H
#define CGAL_FILL_HOLE_POLYHEDRON_H
#include <CGAL/internal/Fair.h>
#include <CGAL/internal/Refine.h>
#include <CGAL/internal/Triangulate.h>

#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <CGAL/internal/Weights.h>
#include <vector>
#include <limits>
#include <set>
#include <map>

// for default parameters
#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>  // for sparse linear system solver
#if defined(CGAL_SUPERLU_ENABLED)
#include <Eigen/SuperLUSupport>
#else
#include <Eigen/SparseLU>
#endif
#endif

namespace CGAL {

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

  void average_length(Polyhedron& poly, Vertex_handle vh, std::map<Vertex_handle, double>& scale_attribute)
  {
    const Point_3& vp = vh->point(); 
    Halfedge_around_vertex_circulator circ(vh->vertex_begin()), done(circ);
    int deg = 0;
    double sum = 0;
    do {
      const Point_3& vq = circ->opposite()->vertex()->point();
      sum += sqrt(CGAL::squared_distance(vp, vq));
      ++deg;
      ++circ;
    } while(circ != done);
    scale_attribute[vh] = sum/deg;
  }

public:
  void operator()(Polyhedron& poly, Halfedge_handle it)
  {
    std::cerr << "fill" << std::endl;
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
    triangulate_functor(it, poly, facets);
    
    if(is_refining) {
      // Refine
      internal::Refine_Polyhedron_3<Polyhedron> refine_functor(alpha);
      refine_functor(scale_attribute, facets, poly);

      std::cerr << "|facets| = " << facets.size() << std::endl;

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

        std::cerr << "|boundary vertices| = " << boundary_vertices.size() << std::endl;
        std::cerr << "|interior vertices| = " << interior_vertices.size() << std::endl;
        // Fair
        internal::Fair_Polyhedron_3<Polyhedron, SparseLinearSolver, WeightCalculator> fair_functor(weight_calculator);
        fair_functor(interior_vertices, poly);
      } // if(is_fairing)
    } // if(is_refining)
  }
};

enum Fairing_weight_type_tag {
  UNIFORM_WEIGHTING,
  SCALE_DEPENDENT_WEIGHTING,
  COTANGENT_WEIGHTING
};

template<class Sparse_linear_system, class Polyhedron>
void fill(Polyhedron& poly, 
  typename Polyhedron::Halfedge_handle it, 
  bool refine = true,
  double density_control_factor = 1.41 /* ~sqrt(2) */,
  bool fair = true,
  Fairing_weight_type_tag weight_tag = SCALE_DEPENDENT_WEIGHTING
  )
{
  typedef CGAL::internal::Uniform_weight_fairing<Polyhedron> Uniform_weight;
  typedef CGAL::internal::Cotangent_weight_with_voronoi_area_fairing<Polyhedron> Cotangent_weight;
  typedef CGAL::internal::Scale_dependent_weight_fairing<Polyhedron> Scale_dependent_weight;

  if(weight_tag == COTANGENT_WEIGHTING) {
    Fill_hole_Polyhedron_3<Polyhedron, Sparse_linear_system, Cotangent_weight>
      (refine, density_control_factor, fair)(poly, it);
  }
  else if(weight_tag == UNIFORM_WEIGHTING) {
    Fill_hole_Polyhedron_3<Polyhedron, Sparse_linear_system, Uniform_weight>
      (refine, density_control_factor, fair)(poly, it);
  }
  else if(weight_tag == SCALE_DEPENDENT_WEIGHTING) {
    Fill_hole_Polyhedron_3<Polyhedron, Sparse_linear_system, Scale_dependent_weight>
      (refine, density_control_factor, fair)(poly, it);
  }
}

template<class Polyhedron>
void fill(Polyhedron& poly, 
  typename Polyhedron::Halfedge_handle it, 
  bool refine = true,
  double density_control_factor = 1.41 /* ~sqrt(2) */,
  bool fair = true,
  Fairing_weight_type_tag weight_tag = SCALE_DEPENDENT_WEIGHTING
  )
{
  typedef
#if defined(CGAL_EIGEN3_ENABLED)
#if defined(CGAL_SUPERLU_ENABLED)
  CGAL::Eigen_solver_traits<Eigen::SuperLU<CGAL::Eigen_sparse_matrix<double>::EigenType> >
#else
  CGAL::Eigen_solver_traits<
    Eigen::SparseLU<
      CGAL::Eigen_sparse_matrix<double, Eigen::ColMajor>::EigenType,
      Eigen::COLAMDOrdering<int> >  >
#endif
#endif
  Sparse_linear_system;

  fill<Sparse_linear_system, Polyhedron>(poly, it, refine, density_control_factor, fair, weight_tag);
}


} // namespace CGAL

#endif