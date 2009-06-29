// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de> 
//
// ============================================================================

// Contains:  
// - Test compute_smallest_nonnegative_root
// - Test compare_smallest_nonnegative_roots


#include <CGAL/basic.h>
#include <CGAL/algebraic_kernel_1_tools.h>

template <class Algebraic_kernel_1>
void test_algebraic_kernel_1_tools(){
  typedef typename Algebraic_kernel_1::Coefficient Coefficient;
  typedef typename Algebraic_kernel_1::Boundary Boundary;
  
  typedef typename Algebraic_kernel_1::Algebraic_real_1 Root;
  typedef typename Algebraic_kernel_1::Polynomial_1 Polynomial_1;
 
  typedef CGAL::Polynomial_traits_d<Polynomial_1> PT_1;
  
  Algebraic_kernel_1 ak; 
  
  Polynomial_1 x = typename PT_1::Shift()(Polynomial_1(1),1);
  Polynomial_1 p1 = x*x-2;
  Polynomial_1 p2 = x*x-3;

  typename Algebraic_kernel_1::Solve_1 solve_1 = ak.solve_1_object();
  
  std::vector<Root> roots1, roots2;
  solve_1(p1,std::back_inserter(roots1));
  solve_1(p2,std::back_inserter(roots2));
  
  assert(roots1.size() == 2);
  assert(roots2.size() == 2);
  
  Root r1 = roots1[1]; // sqrt(2)
  Root r2 = roots2[1]; // sqrt(3)
  
  assert(r1 == *CGAL::compute_smallest_nonnegative_root(ak,p1));
  assert(r2 == *CGAL::compute_smallest_nonnegative_root(ak,p2));
  
  assert(CGAL::compare_smallest_nonnegative_roots(ak,p1,p2) 
      == CGAL::compare(r1,r2));                 
  
}


int main(){
  typedef CGAL::Arithmetic_kernel Arithmetic_kernel;
  typedef Arithmetic_kernel::Integer Coefficient;
  typedef Arithmetic_kernel::Rational Boundary;
  
  typedef CGAL::Algebraic_kernel_1<Coefficient, Boundary> Algebraic_kernel_1;
 
  test_algebraic_kernel_1_tools<Algebraic_kernel_1>();
  return 0;
}



