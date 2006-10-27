// Example: computes point of minimum norm in the intersection of halfspaces 
#include <cassert>
#include <vector>
#include <CGAL/Cartesian_d.h>
#include <CGAL/MP_Float.h>

#include <boost/iterator/counting_iterator.hpp>
#include <CGAL/iterator.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// unary function that maps integer j to iterator constructed from j
template <class Iterator>
struct Jth {
  typedef Iterator result_type;
  result_type operator() (int j) const { return Iterator (j); }
}; 

// unary function that maps a hyperplane h to its negative j-th coordinate,
// for some fixed j
template <class Hyperplane>
class Fixed_coordinate {
public:
  typedef typename CGAL::Kernel_traits<Hyperplane>::Kernel::RT result_type;
  Fixed_coordinate (int index) : j (index) {}
  result_type operator() (const Hyperplane& h) const {return -h[j];} 
private:
  const int j;
};

// unary function that maps an integer i to 1 if i=j and 0 otherwise,
// for some fixed j
template <class RT>
class Unit_vector {
public:
  typedef RT result_type;
  Unit_vector (int index) : j (index) {}
  result_type operator() (int i) { return i==j ? 1 : 0; }
private:
  const int j;
};

// function to solve the QP that computes the minimum-norm point in the
// intersection of positive halfspaces induced by a set of oriented 
// hyperplanes 
template <class HyperplaneIterator, class ET>
CGAL::Quadratic_program_solution<ET> 
solve_minimum_norm_point_qp 
(HyperplaneIterator begin, HyperplaneIterator end, const ET& dummy)
{
  // hyperplane type
  typedef typename HyperplaneIterator::value_type H;
  // input number type
  typedef typename CGAL::Kernel_traits<H>::Kernel::RT RT; 
  // iterator type for j-th column of A (j=0,...,d-1) and for b (j=d);
  // transforms iterator for hyperplane to iterator for j-th coordinates
  typedef CGAL::Join_input_iterator_1
    <HyperplaneIterator, Fixed_coordinate<H> > Fixed_coordinate_iterator;
  // iterator type for A; transforms iterator for indices to iterator
  // for corresponding columns of A
  typedef CGAL::Join_input_iterator_1
    <boost::counting_iterator<int>, Jth<Fixed_coordinate_iterator> > A_it;
  // iterator type for b; 
  typedef Fixed_coordinate_iterator B_it;
  // iterator type for relations (all are "<=")
  typedef CGAL::Const_oneset_iterator<CGAL::Comparison_result> R_it;
  // iterator type for c (all entries are 0)
  typedef CGAL::Const_oneset_iterator<RT> C_it;
  // iterator type for i-th row of D; transforms iterator for indices
  // to iterator for the i-th unit vector
  typedef CGAL::Join_input_iterator_1
    <boost::counting_iterator<int>, Unit_vector<RT> > Unit_vector_iterator;
  // iterator type for D; transforms iterator for indices to iterator
  // for corresponding rows of D
  typedef CGAL::Join_input_iterator_1
    <boost::counting_iterator<int>, Jth<Unit_vector_iterator> > D_it;
  // determine dimension from first hyperplane
  int d = (*begin).dimension();
  // construct program and solve itx
  return CGAL::solve_quadratic_program 
    (CGAL::make_free_quadratic_program_from_iterators 
     (d, end-begin, A_it(0), B_it(d), R_it(CGAL::SMALLER), D_it(), C_it (0)));
}

int main() {
  // 
  return 0;

}
