
template < class Point >
class PD_D_iterator {

  // types needed for iterators
  // ...
  typedef  typename std::vector<Point>::const_iterator    Point_it;
  typedef  PD_D_row_iterator<Point>                       Row_it;
  
  // data members
  Point_it      points;       // iterator referring to C matrix
  int           i;            // row index
  
public:
  PD_D_iterator( Point_it it) : points( it), i( 0) { }
  
  Row_it  operator * ( ) { return Row_it( points, i);}
  
  // other operators required for random access iterators
  // ...
  
};


typedef  CGAL::QP_transform_iterator_1< Matrix::const_iterator>  Vector_iterator;

struct QPSolverTraits {
    typedef  GMP::Double                         ET;
    typedef  Vector_iterator                     A_iterator;
    typedef  CGAL::QP_const_value_iterator<IT>  B_iterator;
    typedef  CGAL::QP_const_value_iterator<IT>  C_iterator;
    typedef  PD_D_iterator<Point>                D_iterator;

    enum Row_type { LESS_EQUAL, EQUAL, GREATER_EQUAL};
    typedef  CGAL::QP_const_value_iterator<Row_type>  Row_type_iterator;
    
    typedef  CGAL::Tag_false  Is_linear;
    typedef  CGAL::Tag_true   Is_symmetric;
    typedef  CGAL::Tag_true   Has_full_row_rank;
    typedef  CGAL::Tag_false  Use_perturbation;
};

typedef CGAL::QP_solver<QPSolverTraits>     Solver;

int main( int argc, char** argv)
{
  Matrix  A( 6);
  Matrix  C( 6);

  // constraint matrix A, column by column
  A[ 0].push_back( 1); A[ 1].push_back( 1); A[ 2].push_back( 1);
  A[ 0].push_back( 0); A[ 1].push_back( 0); A[ 2].push_back( 0);
  A[ 3].push_back( 0); A[ 4].push_back( 0); A[ 5].push_back( 0);
  A[ 3].push_back( 1); A[ 4].push_back( 1); A[ 5].push_back( 1);

  // point matrix C, column by column 
  C[ 0].push_back(  2); C[ 1].push_back(  6); C[ 2].push_back(  5);
  C[ 0].push_back(  2); C[ 1].push_back(  1); C[ 2].push_back(  5);
  
  C[ 3].push_back( -7); C[ 4].push_back( -9); C[ 5].push_back( -6);
  C[ 3].push_back( -4); C[ 4].push_back( -7); C[ 5].push_back( -7);      

  // configure partial filtered pricing strategy with default template
  // parameters
  CGAL::QP_partial_filtered_pricing<QPSolverTraits>  strategy;
  
  // solve qp with (explicit) partial filtered pricing strategy
  Solver       qp( 6, 2, Vector_iterator( A.begin()),
                   QPSolverTraits::B_iterator(1.0),
                   QPSolverTraits::C_iterator(0.0),
                   QPSolverTraits::D_iterator(C.begin()),
                   QPSolverTraits::Row_type_iterator(QPSolverTraits::EQUAL),
                   strategy);
  
  // query solution, if qp is optimal				     
  if (qp.status() == Solver::OPTIMAL) {
    Solver::Variable_value_iterator v_it;
    Solver::Variable_value_iterator e_it = qp.variables_value_end(); 
    for (int i = 0, v_it = qp.variables_value_begin(); v_it != e_it; ++v_it, ++i) {
      std::cout << "x[" << i << "]= " << *v_it << std::endl;
    }
  }			   
  return 0;
}
