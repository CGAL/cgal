#ifndef CGAL_HOLE_FILLING_FAIR_POLYHEDRON_3_H
#define CGAL_HOLE_FILLING_FAIR_POLYHEDRON_3_H

#include <map>
#include <set>
#include <CGAL/assertions.h>
#include <CGAL/internal/Hole_filling/Weights.h>
#include <CGAL/Timer.h>
#include <CGAL/trace.h>

// for default parameters
#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>  // for sparse linear system solver
#endif

namespace CGAL {
/*!
\ingroup PkgHoleFilling
@brief Fairing continuity type
*/
enum Fairing_continuity 
{ 
  FAIRING_C_0 = 0, /**< C0 continuity */
  FAIRING_C_1 = 1,  /**< C1 continuity */
  FAIRING_C_2 = 2   /**< C2 continuity */
};

namespace internal {

struct Fair_default_sparse_linear_solver {
  typedef
#if defined(CGAL_EIGEN3_ENABLED)
  #if defined(CGAL_SUPERLU_ENABLED)
    CGAL::Eigen_solver_traits<Eigen::SuperLU<CGAL::Eigen_sparse_matrix<double>::EigenType> >
  #else
    #if EIGEN_VERSION_AT_LEAST(3,2,0)
    CGAL::Eigen_solver_traits<
      Eigen::SparseLU<
        CGAL::Eigen_sparse_matrix<double>::EigenType,
        Eigen::COLAMDOrdering<int> >  >
    #else
    Fair_default_sparse_linear_solver // dummy type to make it compile
    #endif
  #endif
#else
  Fair_default_sparse_linear_solver // dummy type to make it compile
#endif
  Solver;
};

// [On Linear Variational Surface Deformation Methods-2008]
template<class Polyhedron, class SparseLinearSolver, class WeightCalculator>
class Fair_Polyhedron_3 {
// typedefs
  typedef typename Polyhedron::Traits::Point_3 Point_3;
  typedef typename Polyhedron::Vertex_handle Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;

  typedef SparseLinearSolver Sparse_linear_solver;
// members
  WeightCalculator weight_calculator;
public:
  Fair_Polyhedron_3(WeightCalculator weight_calculator = WeightCalculator())
    : weight_calculator(weight_calculator) { }

private:
  double sum_weight(Vertex_handle v, Polyhedron& polyhedron) {
    double weight = 0;
    Halfedge_around_vertex_circulator circ = v->vertex_begin();
    do {
      weight += weight_calculator.w_ij(circ, polyhedron);
    } while(++circ != v->vertex_begin());
    return weight;
  }

  // recursively computes a row (use depth parameter to compute L, L^2, L^3)
  // Equation 6 in [On Linear Variational Surface Deformation Methods]
  void compute_row(
    Vertex_handle v, 
    Polyhedron& polyhedron,
    std::size_t row_id,                            // which row to insert in [ frees stay left-hand side ]
    typename Sparse_linear_solver::Matrix& matrix, 
    double& x, double& y, double& z,               // constants transfered to right-hand side
    double multiplier,
    const std::map<Vertex_handle, std::size_t>& vertex_id_map,
    unsigned int depth)
  {
    if(depth == 0) {
      typename std::map<Vertex_handle, std::size_t>::const_iterator vertex_id_it = vertex_id_map.find(v);
      if(vertex_id_it != vertex_id_map.end()) {
        matrix.add_coef(row_id, vertex_id_it->second, multiplier);
      }
      else { 
        x += multiplier * -v->point().x(); 
        y += multiplier * -v->point().y(); 
        z += multiplier * -v->point().z(); 
      }
    }
    else {
      double w_i = weight_calculator.w_i(v, polyhedron);

      Halfedge_around_vertex_circulator circ = v->vertex_begin();
      do {
        double w_i_w_ij = w_i * weight_calculator.w_ij(circ, polyhedron) ;

        Vertex_handle nv = circ->opposite()->vertex();
        compute_row(nv, polyhedron, row_id, matrix, x, y, z, -w_i_w_ij*multiplier, vertex_id_map, depth-1);
      } while(++circ != v->vertex_begin());

      double w_i_w_ij_sum = w_i * sum_weight(v, polyhedron);
      compute_row(v, polyhedron, row_id, matrix, x, y, z, w_i_w_ij_sum*multiplier, vertex_id_map, depth-1);
    }
  }

public:
  template<class InputIterator>
  bool fair(Polyhedron& polyhedron, InputIterator vb, InputIterator ve, Fairing_continuity fc)
  {
    int depth = static_cast<int>(fc) + 1;
    if(depth < 1 || depth > 3) {
      CGAL_warning(!"Continuity should be between 0 and 2 inclusively!");
      return false; 
    }

    std::set<Vertex_handle> interior_vertices(vb, ve);
    if(interior_vertices.empty()) { return true; }

    CGAL::Timer timer; timer.start();

    const std::size_t nb_vertices = interior_vertices.size();
    typename Sparse_linear_solver::Vector X(nb_vertices), Bx(nb_vertices);
    typename Sparse_linear_solver::Vector Y(nb_vertices), By(nb_vertices);
    typename Sparse_linear_solver::Vector Z(nb_vertices), Bz(nb_vertices);

    std::map<Vertex_handle, std::size_t> vertex_id_map;
    std::size_t id = 0;
    for(typename std::set<Vertex_handle>::iterator it = interior_vertices.begin(); it != interior_vertices.end(); ++it, ++id) {
      if( !vertex_id_map.insert(std::make_pair(*it, id)).second ) {
        CGAL_warning(!"Duplicate vertex is found!");
        return false;
      }
    }

    typename Sparse_linear_solver::Matrix A(nb_vertices);

    for(typename std::set<Vertex_handle>::iterator vb = interior_vertices.begin(); vb != interior_vertices.end(); ++vb) {
      std::size_t v_id = vertex_id_map[*vb];
      compute_row(*vb, polyhedron, v_id, A, Bx[v_id], By[v_id], Bz[v_id], 1, vertex_id_map, depth);
    }
    CGAL_TRACE_STREAM << "**Timer** System construction: " << timer.time() << std::endl; timer.reset();

    // factorize
    double D;
    Sparse_linear_solver m_solver;
    bool prefactor_ok = m_solver.factor(A, D);
    if(!prefactor_ok) {
      CGAL_warning(!"pre_factor failed!");
      return false;
    }
    CGAL_TRACE_STREAM << "**Timer** System factorization: " << timer.time() << std::endl; timer.reset();

    // solve
    bool is_all_solved = m_solver.linear_solver(Bx, X) && m_solver.linear_solver(By, Y) && m_solver.linear_solver(Bz, Z);
    if(!is_all_solved) {
      CGAL_warning(!"linear_solver failed!"); 
      return false; 
    }
    CGAL_TRACE_STREAM << "**Timer** System solver: " << timer.time() << std::endl; timer.reset();

    
    /* This relative error is to large for cases that the results are not good */ 
    /*
    double rel_err_x = (A.eigen_object()*X - Bx).norm() / Bx.norm();
    double rel_err_y = (A.eigen_object()*Y - By).norm() / By.norm();
    double rel_err_z = (A.eigen_object()*Z - Bz).norm() / Bz.norm();
    CGAL_TRACE_STREAM << "rel error: " << rel_err_x 
                                << " " << rel_err_y
                                << " " << rel_err_z << std::endl;
                                */

    // update 
    id = 0;
    for(typename std::set<Vertex_handle>::iterator it = interior_vertices.begin(); it != interior_vertices.end(); ++it, ++id) {
      (*it)->point() = Point_3(X[id], Y[id], Z[id]);
    }
    return true;
  }
};

}//namespace internal

/*!
\ingroup PkgHoleFilling
@brief Function fairing a region on surface mesh. 
The region denoted by @a vertex_begin and @a vertex_end might contain multiple disconnected components.
Note that the structure is not altered in any way, only positions of the vertices get updated.

Fairing might fail if fixed vertices, which are used as boundary conditions, do not suffice to solve constructed linear system.
The larger @a continuity gets, the more fixed vertices are required.

@tparam SparseLinearSolver a model of SparseLinearAlgebraTraitsWithPreFactor_d. If \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available 
        and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits` is provided as default parameter.
@tparam WeightCalculator a model of FairWeightCalculator and can be omitted to use default Cotangent weights
@tparam Polyhedron a %CGAL polyhedron
@tparam InputIterator iterator over input vertices

@param polyhedron surface mesh to be faired
@param vertex_begin first iterator of the range of vertices
@param vertex_end past-the-end iterator of the range of vertices
@param weight_calculator function object to calculate weights, default to Cotangent weights and can be omitted
@param continuity tangential continuity, default to FAIRING_C_1 and can be omitted

@return true if fairing is successful, otherwise no vertex position is changed

@todo accuracy of solvers are not good, for example when there is no boundary condition pre_factor should fail, but it does not.
*/
template<class SparseLinearSolver, class WeightCalculator, class Polyhedron, class InputIterator>
bool fair(Polyhedron& polyhedron, 
  InputIterator vertex_begin,
  InputIterator vertex_end,
  WeightCalculator weight_calculator,
  Fairing_continuity continuity = FAIRING_C_1
  )
{
  internal::Fair_Polyhedron_3<Polyhedron, SparseLinearSolver, WeightCalculator> fair_functor(weight_calculator);
  return fair_functor.fair(polyhedron, vertex_begin, vertex_end, continuity);
}

//use default SparseLinearSolver
template<class WeightCalculator, class Polyhedron, class InputIterator>
bool fair(Polyhedron& poly, 
  InputIterator vb,
  InputIterator ve,
  WeightCalculator weight_calculator,
  Fairing_continuity continuity = FAIRING_C_1
  )
{
  typedef internal::Fair_default_sparse_linear_solver::Solver Sparse_linear_solver;
  return fair<Sparse_linear_solver, WeightCalculator, Polyhedron, InputIterator>
    (poly, vb, ve, weight_calculator, continuity);
}

//use default WeightCalculator
template<class SparseLinearSolver, class Polyhedron, class InputIterator>
bool fair(Polyhedron& poly, 
  InputIterator vb,
  InputIterator ve,
  Fairing_continuity continuity = FAIRING_C_1
  )
{
  typedef internal::Cotangent_weight_with_voronoi_area_fairing<Polyhedron> Weight_calculator;
  return fair<SparseLinearSolver, Weight_calculator, Polyhedron, InputIterator>
    (poly, vb, ve, Weight_calculator(), continuity);
}

//use default SparseLinearSolver and WeightCalculator
template<class Polyhedron, class InputIterator>
bool fair(Polyhedron& poly, 
  InputIterator vb,
  InputIterator ve,
  Fairing_continuity continuity = FAIRING_C_1
  )
{
  typedef internal::Fair_default_sparse_linear_solver::Solver Sparse_linear_solver;
  return fair<Sparse_linear_solver, Polyhedron, InputIterator>(poly, vb, ve, continuity);
}

}//namespace CGAL
#endif //CGAL_HOLE_FILLING_FAIR_POLYHEDRON_3_H
