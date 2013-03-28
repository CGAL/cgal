// example: function to check whether a point is in the convex 
// hull of other points; this version uses a maker
#include <boost/iterator/transform_iterator.hpp>
#include <CGAL/Kernel_traits.h>
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

// function to test whether point is in the convex hull of other points;
// the type ET is an exact type used for the computations
template <class Point_d, class RandomAccessIterator, class ET>
CGAL::Quadratic_program_solution<ET>
solve_convex_hull_containment_lp (const Point_d& p,
				  RandomAccessIterator begin,
				  RandomAccessIterator end, const ET& dummy)
{
  // construct program and solve it
  return CGAL::solve_nonnegative_linear_program
    (CGAL::make_nonnegative_linear_program_from_iterators
     (static_cast<int>(end-begin),                                         // n
      p.dimension()+1,                                                     // m
      boost::transform_iterator
      <Homogeneous_begin<Point_d>, RandomAccessIterator>(begin),           // A
      typename Point_d::Homogeneous_const_iterator (p.homogeneous_begin()),// b
      CGAL::Const_oneset_iterator<CGAL::Comparison_result>(CGAL::EQUAL),   // ~
      CGAL::Const_oneset_iterator
      <typename CGAL::Kernel_traits<Point_d>::Kernel::RT> (0)), dummy);    // c
}
