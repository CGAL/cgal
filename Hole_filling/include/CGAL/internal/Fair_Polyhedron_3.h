#ifndef CGAL_HOLE_FILLING_FAIR_POLYHEDRON_3_H
#define CGAL_HOLE_FILLING_FAIR_POLYHEDRON_3_H

#include <map>
#include <set>
#include <CGAL/assertions.h>
#include <CGAL/internal/Weights.h>
#include <CGAL/Timer.h>
#include <CGAL/trace.h>

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
namespace internal {

struct Fair_default_sparse_linear_solver {
  typedef
#if defined(CGAL_EIGEN3_ENABLED)
  #if defined(CGAL_SUPERLU_ENABLED)
    CGAL::Eigen_solver_traits<Eigen::SuperLU<CGAL::Eigen_sparse_matrix<double>::EigenType> >
  #else
    CGAL::Eigen_solver_traits<
      Eigen::SparseLU<
        CGAL::Eigen_sparse_matrix<double>::EigenType,
        Eigen::COLAMDOrdering<int> >  >
  #endif
#else
  Fair_default_sparse_linear_solver // dummy type to make it compile
#endif
  Solver;
};

template<class SparseMatrix>
void write_matrix_sparse(std::ostream &stream, const SparseMatrix& m)
{
  stream << m.rows() << " " << m.cols() << std::endl;

  for(typename SparseMatrix::Index r = 0; r < m.rows(); ++r) 
  {
    for(typename SparseMatrix::Index c = 0; c < m.cols(); ++c) 
    {
      //if(m.coeff(r,c) != 0)
      stream << m.coeff(r,c) << " ";
    }
    stream << std::endl;
  }
}

template<class Matrix>
void write_matrix(std::ostream &stream, const Matrix& m)
{
  stream << m.rows() << " " << m.cols() << std::endl;

  for(typename Matrix::Index r = 0; r < m.rows(); ++r) 
  {
    for(typename Matrix::Index c = 0; c < m.cols(); ++c) 
    {
      stream << m(r,c) << " ";
    }
    stream << std::endl;
  }
}

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

  void one_ring(Vertex_handle v, 
    Polyhedron& polyhedron, 
    std::vector<std::pair<std::size_t, double> >& row, 
    double& x, double& y, double& z, 
    const std::map<Vertex_handle, std::size_t>& vertex_id_map) 
  {
    typedef typename std::map<Vertex_handle, std::size_t>::const_iterator Vertex_id_it;
    row.clear();
    x = y = z = 0.0;

    double w_v = weight_calculator.w_i(v, polyhedron);
    double diagonal_v = 0.0;
    Halfedge_around_vertex_circulator circ = v->vertex_begin();
    do {
      Vertex_handle nv = circ->opposite()->vertex();          
      double nv_weight = weight_calculator.w_ij(circ, polyhedron);
      double nv_total_weight = w_v * nv_weight ;
      diagonal_v += -nv_total_weight;

      Vertex_id_it vertex_id_it = vertex_id_map.find(nv);
      if(vertex_id_it != vertex_id_map.end()) {
        row.push_back(std::make_pair(vertex_id_it->second, nv_total_weight));
      }
      else {
        x += - nv_total_weight * nv->point().x();
        y += - nv_total_weight * nv->point().y();
        z += - nv_total_weight * nv->point().z();
      }
    } while(++circ != v->vertex_begin());

    Vertex_id_it vertex_id_it = vertex_id_map.find(v);
    if(vertex_id_it != vertex_id_map.end()) {
      row.push_back(std::make_pair(vertex_id_it->second, diagonal_v));
    }
    else {
      x += -diagonal_v * v->point().x();
      y += -diagonal_v * v->point().y();
      z += -diagonal_v * v->point().z();
    }
  }

public:
  template<class InputIterator>
  bool fair(Polyhedron& polyhedron, InputIterator vb, InputIterator ve)
  {
    std::set<Vertex_handle> interior_vertices(vb, ve);
    if(interior_vertices.empty()) { return false; }

    CGAL::Timer timer; timer.start();

    const std::size_t nb_vertices = interior_vertices.size();
    typename Sparse_linear_solver::Vector X(nb_vertices), Bx(nb_vertices);
    typename Sparse_linear_solver::Vector Y(nb_vertices), By(nb_vertices);
    typename Sparse_linear_solver::Vector Z(nb_vertices), Bz(nb_vertices);

    std::map<Vertex_handle, std::size_t> vertex_id_map;
    std::size_t id = 0;
    for(typename std::set<Vertex_handle>::iterator it = interior_vertices.begin(); it != interior_vertices.end(); ++it, ++id) {
      vertex_id_map[*it] = id;      
    }

    typename Sparse_linear_solver::Matrix A(nb_vertices);
    std::vector<std::pair<std::size_t, double> > row;

    for(typename std::set<Vertex_handle>::iterator vb = interior_vertices.begin(); vb != interior_vertices.end(); ++vb) {
      std::size_t v_id = vertex_id_map[*vb];
      double w_v = weight_calculator.w_i(*vb, polyhedron);
      double wn_sum = sum_weight(*vb, polyhedron);
      double x, y, z;
      x = y = z = 0.0;
      double diagonal_v = 0.0;

      double nx, ny, nz;
      nx = ny = nz = 0;
      one_ring(*vb, polyhedron, row, nx, ny, nz, vertex_id_map);
      for(std::vector<std::pair<std::size_t, double> >::iterator it = row.begin(); it != row.end(); ++it) {
        A.add_coef(v_id, it->first, - w_v * wn_sum * it->second);
      }
      x += - w_v * wn_sum * nx; 
      y += - w_v * wn_sum * ny; 
      z += - w_v * wn_sum * nz;

      Halfedge_around_vertex_circulator circ = (*vb)->vertex_begin();
      do {
        Vertex_handle nv = circ->opposite()->vertex();
        double nv_w = weight_calculator.w_ij(circ, polyhedron);

        one_ring(nv, polyhedron, row, nx, ny, nz, vertex_id_map);
        for(std::vector<std::pair<std::size_t, double> >::iterator it = row.begin(); it != row.end(); ++it) {
          A.add_coef(v_id, it->first, w_v * nv_w * it->second);
        }
        x += w_v * nv_w * nx; 
        y += w_v * nv_w * ny; 
        z += w_v * nv_w * nz;
      } while(++circ != (*vb)->vertex_begin());

      Bx[v_id] = x; By[v_id] = y; Bz[v_id] = z;
    }

    CGAL_TRACE_STREAM << "**Timer** System construction: " << timer.time() << std::endl; timer.reset();

    //std::ofstream out_A("D:/A.txt");
    //std::ofstream out_b("D:/b.txt");
    //write_matrix_sparse(out_A, A.eigen_object());
    //write_matrix(out_b, Bx);

    // factorize
    double D;
    Sparse_linear_solver m_solver;
    bool prefactor_ok = m_solver.pre_factor(A, D);
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
    // update 
    id = 0;
    for(typename std::set<Vertex_handle>::iterator it = interior_vertices.begin(); it != interior_vertices.end(); ++it, ++id) {
      (*it)->point() = Point_3(X[id], Y[id], Z[id]);
    }
    return true;
  }
};

}//namespace internal

/** 
 * @brief Function fairing a region on surface mesh.
 *
 * @tparam SparseLinearSolver a model of SparseLinearAlgebraTraitsWithPreFactor_d and can be omitted if Eigen defined...(give exact models etc)
 * @tparam WeightCalculator a model of "weight model" and can be omitted to use default Cotangent weights
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam InputIterator iterator over input vertices
 *
 * @param polyhedron surface mesh to be faired
 * @param vertex_begin first iterator of the range of vertices
 * @param vertex_end past-the-end iterator of the range of vertices
 * @param weight_calculator function object to calculate weights
 */
template<class SparseLinearSolver, class WeightCalculator, class Polyhedron, class InputIterator>
bool fair(Polyhedron& polyhedron, 
  InputIterator vertex_begin,
  InputIterator vertex_end,
  WeightCalculator weight_calculator
  )
{
  internal::Fair_Polyhedron_3<Polyhedron, SparseLinearSolver, WeightCalculator> fair_functor(weight_calculator);
  return fair_functor.fair(polyhedron, vertex_begin, vertex_end);
}

//use default SparseLinearSolver
template<class WeightCalculator, class Polyhedron, class InputIterator>
bool fair(Polyhedron& poly, 
  InputIterator vb,
  InputIterator ve,
  WeightCalculator weight_calculator
  )
{
  typedef CGAL::internal::Fair_default_sparse_linear_solver::Solver Sparse_linear_solver;
  return fair<Sparse_linear_solver, WeightCalculator, Polyhedron, InputIterator>(poly, vb, ve, weight_calculator);
}

//use default WeightCalculator
template<class SparseLinearSolver, class Polyhedron, class InputIterator>
bool fair(Polyhedron& poly, 
  InputIterator vb,
  InputIterator ve
  )
{
  typedef CGAL::Cotangent_weight_with_voronoi_area_fairing<Polyhedron> Weight_calculator;
  return fair<SparseLinearSolver, Weight_calculator, Polyhedron, InputIterator>(poly, vb, ve, Weight_calculator());
}

//use default SparseLinearSolver and WeightCalculator
template<class Polyhedron, class InputIterator>
bool fair(Polyhedron& poly, 
  InputIterator vb,
  InputIterator ve
  )
{
  typedef CGAL::internal::Fair_default_sparse_linear_solver::Solver Sparse_linear_solver;
  return fair<Sparse_linear_solver, Polyhedron, InputIterator>(poly, vb, ve);
}

}//namespace CGAL
#endif //CGAL_HOLE_FILLING_FAIR_POLYHEDRON_3_H