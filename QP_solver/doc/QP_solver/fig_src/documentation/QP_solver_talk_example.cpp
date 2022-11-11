
#include <CGAL/QP_solver.h>
#include <CGAL/QP_full_exact_pricing.h>
#include <CGAL/QP_partial_exact_pricing.h>
#include <CGAL/QP_full_filtered_pricing.h>
#include <CGAL/QP_partial_filtered_pricing.h>
#include <CGAL/_QP_solver/Double.h>
#include <iostream>

typedef  std::vector<double>  Vector;
typedef  std::vector<Vector>  Matrix;

typedef CGAL::QP_transform_iterator_1<Matrix::const_iterator>
                                                 Vector_iterator;
typedef Vector::const_iterator                   Entry_iterator;

struct QPSolverTraits {
    typedef  GMP::Double      ET;
    typedef  Vector_iterator  A_iterator;
    typedef  Entry_iterator   B_iterator;
    typedef  Entry_iterator   C_iterator;
    typedef  Vector_iterator  D_iterator;

    enum Row_type { LESS_EQUAL, EQUAL, GREATER_EQUAL};
    typedef  CGAL::QP_const_value_iterator<Row_type>  Row_type_iterator;

    typedef  CGAL::Tag_false  Is_linear;
    typedef  CGAL::Tag_false  Is_symmetric;
    typedef  CGAL::Tag_true   Has_full_row_rank;
    typedef  CGAL::Tag_false  Use_perturbation;
};

typedef CGAL::QP_solver<QPSolverTraits>   Solver;





int main( int argc, char** argv)
{
  Matrix  A( 3);
  Vector  b,c;
  Matrix  D( 3);

  // constraint matrix A, column by column
  A[ 0].push_back( -4); A[ 1].push_back(  2); A[ 2].push_back(  0);
  A[ 0].push_back(  1); A[ 1].push_back(  1); A[ 2].push_back(  1);

  // column vector b
  b.push_back( -8); b.push_back(  2);

  // linear part of objective function, vector c
  c.push_back(  0); c.push_back(  5); c.push_back(  0);

  // quadratic part of objective function, matrix D, row by row
  D[ 0].push_back(  64); D[ 0].push_back( -16); D[ 0].push_back(   0);
  D[ 1].push_back( -16); D[ 1].push_back(   4); D[ 1].push_back(   2);
  D[ 2].push_back(   0); D[ 2].push_back(  -2); D[ 2].push_back(   0);

  // solve qp with pricing strategy defaulting to full exact pricing
  // strategy
  Solver  qp(3, 2, Vector_iterator(A.begin()), b.begin(), c.begin(),
             Vector_iterator(D.begin()), QPSolverTraits::Row_type_iterator
             (QPSolverTraits::GREATER_EQUAL));

  // query solution, if qp is optimal
  if (qp.status() == Solver::OPTIMAL) {
    Solver::Variable_value_iterator v_it;
    Solver::Variable_value_iterator e_it = qp.variables_value_end();
    for(int i=0, v_it=qp.variables_value_begin(); v_it!=e_it; ++v_it, ++i){
      std::cout << "x[" << i << "]= " << *v_it << std::endl;
    }
  }
  return 0;
}
