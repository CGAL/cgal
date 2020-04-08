// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
//
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#include <CGAL/Algebraic_kernel_d/flags.h>

#include <vector>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d_1.h>

#include <CGAL/Algebraic_kernel_d/algebraic_curve_kernel_2_tools.h>

template<typename ArithmeticKernel>
void test_routine() {
  typedef ArithmeticKernel Arithmetic_kernel;
  typedef typename Arithmetic_kernel::Integer Integer;
  typedef typename Arithmetic_kernel::Rational Rational;

  typedef typename CGAL::Polynomial_type_generator<Integer,1>::Type Poly_int1;

  typedef CGAL::Algebraic_kernel_d_1<Integer> Algebraic_kernel_d_1;

  Algebraic_kernel_d_1 ak_1;

  typedef typename Algebraic_kernel_d_1::Solve_1 Real_roots;

  typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Algebraic_real;


  {
    Poly_int1 p(1,-3,1,0,-1,0,0,1);
    Real_roots rr;
    std::vector<Algebraic_real> v1;

    rr(p,std::back_inserter(v1),true);

    std::ptrdiff_t n = std::distance(v1.begin(),v1.end());

    assert(n==3);

    std::vector<Rational> inter;
    CGAL::internal::find_intermediate_values
      (&ak_1,v1.begin(),v1.end(),std::back_inserter(inter));
    int m = (int)inter.size();
    assert(m==n+1);

    for(int i=0;i<m-1;i++) {
      assert(v1[i].compare(inter[i])==CGAL::LARGER);
    }
    for(int i=1;i<m;i++) {
      assert(v1[i-1].compare(inter[i])==CGAL::SMALLER);
    }
  }
  {
    Poly_int1 p(0,2,-6,1,2);
    Real_roots rr2;
    std::vector<Algebraic_real> v2;

    rr2(p,std::back_inserter(v2),false);

    std::ptrdiff_t n = std::distance(v2.begin(), v2.end());

    assert(n==4);

    std::vector<Rational> inter;
    CGAL::internal::find_intermediate_values
      (&ak_1,v2.begin(),v2.end(),std::back_inserter(inter));
    int m = (int)inter.size();
    assert(m==n+1);

    for(int i=0;i<m-1;i++) {
      assert(v2[i].compare(inter[i])==CGAL::LARGER);
    }
    for(int i=1;i<m;i++) {
      assert(v2[i-1].compare(inter[i])==CGAL::SMALLER);
    }
  }
  {
    Poly_int1 p(1);
    Real_roots rr1;
    std::vector<Algebraic_real> v1;

    rr1(p,std::back_inserter(v1),false);

    std::ptrdiff_t n = std::distance(v1.begin(), v1.end());

    assert(n==0);

    std::vector<Rational> inter;
    CGAL::internal::find_intermediate_values
      (&ak_1,v1.begin(),v1.end(),std::back_inserter(inter));
    int m = (int)inter.size();
    assert(m==n+1);

    for(int i=0;i<m-1;i++) {
      assert(v1[i].compare(inter[i])==CGAL::LARGER);
    }
    for(int i=1;i<m;i++) {
      assert(v1[i-1].compare(inter[i])==CGAL::SMALLER);
    }
  }
  {
    Poly_int1 p(2,3);
    Real_roots rr2;
    std::vector<Algebraic_real> v2;

    rr2(p,std::back_inserter(v2),false);

    std::ptrdiff_t n = std::distance(v2.begin(), v2.end());

    assert(n==1);

    std::vector<Rational> inter;
    CGAL::internal::find_intermediate_values
      (&ak_1,v2.begin(),v2.end(),std::back_inserter(inter));
    int m = (int)inter.size();
    assert(m==n+1);

    for(int i=0;i<m-1;i++) {
      assert(v2[i].compare(inter[i])==CGAL::LARGER);
    }
    for(int i=1;i<m;i++) {
      assert(v2[i-1].compare(inter[i])==CGAL::SMALLER);
    }
  }
  /* TODO: Move this test into _test_algebraic_kernel_1.h ?
  {
    Algebraic_real a(Integer(2));

    assert(a.is_rational());

    std::stringstream ss("P[39(0,-54465841152)(1,-63020685312)(2,-175239189504)(3,2165206298112)(4,-3645852374016)(5,760508507520)(6,3570449785920)(7,-3732467341536)(8,982269209160)(9,343483161984)(10,-518868252576)(11,1222906488192)(12,-1299164325408)(13,163485500832)(14,760313409786)(15,-718962175440)(16,182270429568)(17,201251093040)(18,-233120831418)(19,114457283784)(20,-13739761548)(21,-28504417464)(22,23634422568)(23,-8334674856)(24,-171675198)(25,2607514920)(26,-1479476808)(27,288652968)(28,29925054)(29,-111025584)(30,75042306)(31,-18475488)(32,5124366)(33,-2154672)(34,205506)(35,-28440)(36,2214)(37,4320)(38,1260)(39,288)]");

    Poly_int1 f;
    ss >> f;

    ::CGAL::set_pretty_mode(std::cout);
    //std::cout << f << " " << a << std::endl;

    assert(CGAL::internal::is_root_of(a,f));

  }
  */
  return;
}

int main() {

#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
  std::cerr << "test LEDA " << std::endl;
  test_routine<CGAL::LEDA_arithmetic_kernel> ();
  std::cerr << "done " << std::endl;
#else
  std::cerr << "test LEDA skipped" << std::endl;
#endif
#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
  std::cerr << "test CORE " << std::endl;
  test_routine<CGAL::CORE_arithmetic_kernel> ();
  std::cerr << "done " << std::endl;
#else
  std::cerr << "test CORE skipped" << std::endl;
#endif
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
  std::cerr << "test GMP " << std::endl;
  test_routine<CGAL::GMP_arithmetic_kernel> ();
  std::cerr << "done " << std::endl;
#else
  std::cerr << "test GMP skipped" << std::endl;
#endif



  return 0;
}
