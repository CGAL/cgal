// example: function to check whether a point is in the convex
// hull of other points
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// unary function to get homogeneous begin-iterator of point
template <class Point_d>
struct Homogeneous_begin  {
  typedef typename Point_d::Homogeneous_const_iterator result_type;
  result_type operator() (const Point_d& p) const {
    return p.homogeneous_begin();
  }
};

// function to solve the LP that tests whether a point is in the
// convex hull of other points; the type ET is an exact type used
// for the internal computations
template <class Point_d, class RandomAccessIterator, class ET>
CGAL::Quadratic_program_solution<ET>
solve_convex_hull_containment_lp (const Point_d& p,
                                  RandomAccessIterator begin,
                                  RandomAccessIterator end, const ET& dummy)
{
  // Constraint matrix type: A[j][i] is the i-th homogeneous coordinate of p_j
  typedef boost::transform_iterator
    <Homogeneous_begin<Point_d>, RandomAccessIterator> A_it;
  // Right-hand side type: b[i] is the i-th homogeneous coordinate of p
  typedef typename Point_d::Homogeneous_const_iterator B_it;
  // Relation type ("=")
  typedef CGAL::Const_oneset_iterator<CGAL::Comparison_result> R_it;
  // input number type
  typedef typename CGAL::Kernel_traits<Point_d>::Kernel::RT RT;
  // Linear objective function type (c=0: we only test feasibility)
  typedef CGAL::Const_oneset_iterator<RT> C_it;
  // the nonnegative linear program type
  typedef
    CGAL::Nonnegative_linear_program_from_iterators<A_it, B_it, R_it, C_it>
    Program;

  // ok, we are prepared now: construct program and solve it
  Program lp (static_cast<int>(end-begin), // number of variables
              p.dimension()+1,             // number of constraints
              A_it (begin), B_it (p.homogeneous_begin()),
              R_it (CGAL::EQUAL), C_it (0));
  return CGAL::solve_nonnegative_linear_program (lp, dummy);
}
