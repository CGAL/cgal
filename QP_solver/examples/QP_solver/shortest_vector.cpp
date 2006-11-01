// Example: computes point of minimum norm in the intersection of halfspaces 
#include <cassert>
#include <vector>
#include <functional>
#include <CGAL/Cartesian_d.h>
#include <CGAL/MP_Float.h>

#include <boost/iterator/counting_iterator.hpp>
#include <CGAL/iterator.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// unary function that maps a hyperplane h to its (negative) j-th coordinate,
// for some fixed j (used to construct iterators for columns of A, and for b) 
template <class Hyperplane>
class Fixed_coordinate {
public:
  typedef typename CGAL::Kernel_traits<Hyperplane>::Kernel::RT result_type;
  Fixed_coordinate (int index, bool negative = false) 
    : j (index), n (negative) {}
  result_type operator() (const Hyperplane& h) const {
    return n ? -h[j] : h[j];} 
private:
  const int j;
  const bool n;
};

// unary function that maps an integer i to 1 if i=j and 0 otherwise,
// for some fixed j (used to construct iterators for rows of D) 
template <class RT>
class Unit_vector {
public:
  typedef RT result_type;
  Unit_vector (int index) : j (index) {}
  result_type operator() (int i) { return i==j ? 1 : 0; }
private:
  const int j;
};

// unary function that maps integer j to iterator constructed from j
// (used to construct iterators for A and D) 
template <class Iterator>
struct Jth {
  typedef Iterator result_type;
  result_type operator() (int j) const { return Iterator (j); }
}; 

// unary function that maps number z to -2z (used to construct iterator for c)
template <class RT>
struct Times_minus_two {
  typedef RT result_type;
  result_type operator() (RT z) const { return -2*z; }
}; 

// function to solve the QP that computes the point closest to a given point
// in the intersection of positive halfspaces induced by a set of oriented 
// hyperplanes; the hyperplane iterator must be random access 
template <class Point_d, class HyperplaneIterator, class ET>
CGAL::Quadratic_program_solution<ET> 
solve_minimum_norm_point_qp 
(const Point_d& p,
 HyperplaneIterator begin, HyperplaneIterator end, const ET& dummy)
{
  // hyperplane type
  typedef typename HyperplaneIterator::value_type H;

  // input number type
  typedef typename CGAL::Kernel_traits<H>::Kernel::RT RT; 

  // iterator type for j-th column of A (j=0,...,d-1) and for b (j=d)
  // transforms iterator for hyperplane to iterator for j-th coordinates
  typedef CGAL::Join_input_iterator_1
    <HyperplaneIterator, Fixed_coordinate<H> > 
    Fixed_coordinate_iterator;

  // iterator type for A; transforms iterator for indices to iterator
  // for corresponding columns of A
  typedef CGAL::Join_input_iterator_1
    <boost::counting_iterator<int>, Jth<Fixed_coordinate_iterator> > A_it;

  // iterator type for relations (all are ">=")
  typedef CGAL::Const_oneset_iterator<CGAL::Comparison_result> R_it;

  // iterator type for i-th row of D; transforms iterator for indices
  // to iterator for the i-th unit vector
  typedef CGAL::Join_input_iterator_1
    <boost::counting_iterator<int>, Unit_vector<RT> > Unit_vector_iterator;

  // iterator type for D; transforms iterator for indices to iterator
  // for corresponding rows of D
  typedef CGAL::Join_input_iterator_1
    <boost::counting_iterator<int>, Jth<Unit_vector_iterator> > D_it;

  // iterator type for c 
  typedef CGAL::Join_input_iterator_1
    <typename Point_d::Cartesian_const_iterator, Times_minus_two<RT> > C_it;

  // now we're set: determine dimension from p, construct program and solve it 
  int d = p.dimension();
  return CGAL::solve_quadratic_program 
    (CGAL::make_free_quadratic_program_from_iterators 
     (d, // n
      end-begin, // m
      A_it(boost::make_counting_iterator<int>(0)), // A,
      Fixed_coordinate_iterator(begin, Fixed_coordinate<H>(d, true)), // b
      R_it(CGAL::LARGER), // ~,
      D_it(boost::make_counting_iterator<int>(0)), // D
      C_it(p.cartesian_begin()) // c
      ), dummy);
}

typedef CGAL::Cartesian_d<double> Kernel_d;
typedef Kernel_d::Point_d Point_d;
typedef Kernel_d::Hyperplane_d Hyperplane_d;

int main() {
  // build the unit square from 4 halfspaces
  std::vector<Hyperplane_d> square;
  square.push_back (Hyperplane_d ( 1,  0,  0));  //  x >= 0
  square.push_back (Hyperplane_d (-1,  0, -1));  //  x <= 1
  square.push_back (Hyperplane_d ( 0,  1,  0));  //  y >= 0
  square.push_back (Hyperplane_d ( 0, -1, -1));  //  y <= 1	
  
  Point_d p(0, 0);

  solve_minimum_norm_point_qp (p, square.begin(), square.end(), 
			       CGAL::MP_Float());
  return 0;

}
