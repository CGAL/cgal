// Example: check whether a point is in the convex hull of other points
#include <cassert>
#include <vector>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/MP_Float.h>
#include <CGAL/iterator.h>
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
bool is_in_convex_hull
(const Point_d& p, 
 RandomAccessIterator begin, RandomAccessIterator end, const ET& dummy) 
{
  // Constraint matrix type: A[j][i] is the i-th homogeneous coordinate of p_j
  typedef CGAL::Join_input_iterator_1
    <RandomAccessIterator, Homogeneous_begin<Point_d> > A_it; 
  // Right-hand side type: b[i] is the i-th homogeneous coordinate of p
  typedef typename Point_d::Homogeneous_const_iterator B_it;
  // Relation type ("=")
  typedef CGAL::Const_oneset_iterator<CGAL::Comparison_result> R_it;
  // input number type
  typedef typename std::iterator_traits<B_it>::value_type RT;
  // Linear objective function type (c=0: we only test feasibility)
  typedef CGAL::Const_oneset_iterator<RT> C_it;
  // the nonnegative linear program type
  typedef
    CGAL::Nonnegative_linear_program_from_iterators<A_it, B_it, R_it, C_it>
    Program;
  // the solution type
  typedef CGAL::Quadratic_program_solution<ET> Solution;

  // ok, we are prepared now: construct program and solve it
  Program lp (end-begin,        // number of variables 
	      p.dimension()+1,  // number of constraints
	      A_it (begin), B_it (p.homogeneous_begin()), 
	      R_it (CGAL::EQUAL), C_it (0));
  Solution s = CGAL::solve_nonnegative_linear_program (lp, ET(0)); 
  
  // p is in the convex hull if and only if the program is feasible
  return (s.status() != CGAL::QP_INFEASIBLE);
}

// now use the above, with MP_Float as the exact type
typedef CGAL::Homogeneous_d<double> Kernel_d;
typedef Kernel_d::Point_d Point_d;

int main()
{
  std::vector<Point_d> points;
  // conex hull: simplex spanned by {(0,0), (10,0), (0,10)} 
  points.push_back (Point_d ( 0.0,  0.0));   
  points.push_back (Point_d (10.0,  0.0));  		    
  points.push_back (Point_d ( 0.0, 10.0));  
  for (int i=0; i<=10; ++i) 
    for (int j=0; j<=10; ++j) {
      // (i,j) is in the simplex iff i+j <= 10 
      bool contained = is_in_convex_hull 
	(Point_d (i, j), points.begin(), points.end(), CGAL::MP_Float());
      assert (contained == (i+j<=10));
    }
	
  return 0;
}
