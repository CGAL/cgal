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

#include <CGAL/Algebraic_curve_kernel_2/flags.h>

#include <CGAL/basic.h>

#include <vector>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_1.h>

#include <CGAL/Algebraic_curve_kernel_2/alg_real_utils.h>

template<typename ArithmeticKernel>
void test_routine() {
  typedef ArithmeticKernel Arithmetic_kernel;
  typedef typename Arithmetic_kernel::Integer Integer;
  typedef typename Arithmetic_kernel::Rational Rational;
  
  typedef typename CGAL::Polynomial_type_generator<Integer,1>::Type Poly_int1;

  typedef CGAL::Algebraic_kernel_1<Integer> Algebraic_kernel_1;
  
  typedef typename Algebraic_kernel_1::Solve_1 Real_roots;

  typedef typename Algebraic_kernel_1::Algebraic_real_1 Algebraic_real;
  
  {
    Poly_int1 p(5,-4,3,-2,-1,1);
    Algebraic_real ar(p,-3,3);
    assert(CGAL::CGALi::estimate_sign_at(ar,p,1000)==CGAL::ZERO);
    Poly_int1 q(1,0,0,-1,1,0,0,1);
    Algebraic_real ar2(q,-2,0);
    assert((CGAL::CGALi::estimate_sign_at(ar2,p))
	       ==CGAL::POSITIVE);
    Poly_int1 r(500001,-400000,300000,-200000,-100000,100000);
    Algebraic_real ar3(r,-5,5);
    assert((CGAL::CGALi::estimate_sign_at(ar3,p))==CGAL::NEGATIVE);
  }
  
  {
    Poly_int1 p(1,-3,1,0,-1,0,0,1);
    Real_roots rr;
    std::vector<Algebraic_real> v1;
    
    rr(p,std::back_inserter(v1),true);
    
    int n = std::distance(v1.begin(),v1.end());

    assert(n==3);
    
    std::vector<Rational> inter;
    CGAL::CGALi::find_intermediate_values(v1.begin(),
                                          v1.end(),
                                          std::back_inserter(inter));
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
    
    rr2(p,std::back_inserter(v2));
    
    int n = std::distance(v2.begin(), v2.end());

    assert(n==4);
    
    std::vector<Rational> inter;
    CGAL::CGALi::find_intermediate_values(v2.begin(),
                                          v2.end(),
                                          std::back_inserter(inter));
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
    
    rr1(p,std::back_inserter(v1));
    
    int n = std::distance(v1.begin(), v1.end());

    assert(n==0);
    
    std::vector<Rational> inter;
    CGAL::CGALi::find_intermediate_values(v1.begin(),
                                          v1.end(),
                                          std::back_inserter(inter));
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
    
    rr2(p,std::back_inserter(v2));
    
    int n = std::distance(v2.begin(), v2.end());

    assert(n==1);
    
    std::vector<Rational> inter;
    CGAL::CGALi::find_intermediate_values(v2.begin(),
                                          v2.end(),
                                          std::back_inserter(inter));
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
    Algebraic_real a(Integer(2));

    assert(a.is_rational());

    std::stringstream ss("P[39(0,-54465841152)(1,-63020685312)(2,-175239189504)(3,2165206298112)(4,-3645852374016)(5,760508507520)(6,3570449785920)(7,-3732467341536)(8,982269209160)(9,343483161984)(10,-518868252576)(11,1222906488192)(12,-1299164325408)(13,163485500832)(14,760313409786)(15,-718962175440)(16,182270429568)(17,201251093040)(18,-233120831418)(19,114457283784)(20,-13739761548)(21,-28504417464)(22,23634422568)(23,-8334674856)(24,-171675198)(25,2607514920)(26,-1479476808)(27,288652968)(28,29925054)(29,-111025584)(30,75042306)(31,-18475488)(32,5124366)(33,-2154672)(34,205506)(35,-28440)(36,2214)(37,4320)(38,1260)(39,288)]");
    
    Poly_int1 f;
    ss >> f; 
    
    ::CGAL::set_pretty_mode(std::cout);
    //std::cout << f << " " << a << std::endl;
    
    assert(CGAL::CGALi::is_root_of(a,f));
    
  }
  {
    Poly_int1 p1(-2,0,1);
    Poly_int1 p2(-3,0,1);
    Algebraic_real sqrt_2(p1,Rational(0),Rational(2));

    assert( CGAL::CGALi::is_root_of(sqrt_2,Poly_int1(0)));
    assert(!CGAL::CGALi::is_root_of(sqrt_2,Poly_int1(1)));
    assert( CGAL::CGALi::is_root_of(sqrt_2,p1));
    assert(!CGAL::CGALi::is_root_of(sqrt_2,p2));  
  }
  
  return;
}

int main(int argc,char** argv) {
  
#ifndef CGAL_USE_LEDA
#ifndef LiS_HAVE_CORE
  std::cerr << "This tests requires LEDA and/or CORE" << std::endl;
  return 1;
#endif
#endif
#ifdef CGAL_USE_LEDA
  // LEDA TEST
  test_routine<CGAL::LEDA_arithmetic_kernel> ();
#else
  std::cerr << "LEDA tests skipped" << std::endl;
#endif
#ifdef CGAL_USE_CORE
  // CORE TEST
  test_routine<CGAL::CORE_arithmetic_kernel> ();
#else
  std::cerr << "CORE tests skipped" << std::endl;
#endif
  return 0;
  
}
