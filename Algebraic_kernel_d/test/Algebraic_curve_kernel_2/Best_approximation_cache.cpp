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

#include <string>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Algebraic_curve_kernel_2/Bitstream_descartes_at_x/Best_approximation_cache.h>

template<typename ArithmeticKernel>
void test_routine() {
  typedef typename ArithmeticKernel::Integer Integer;

  typedef std::string String;

  {
      CGAL::CGALi::Best_approximation_cache<Integer,Integer> approx_mem;

    for(int i=0;i<100;i++) {
      assert(! approx_mem.approximation_exists(Integer(i)));
    }
    int k=1;
    while(k<1000) {
      approx_mem.update_approximation(Integer(k*k),2*k,Integer(k));
      k+=3;
    }
    assert(! approx_mem.approximation_exists(Integer(360)));
    assert(approx_mem.approximation_exists(Integer(361)));
    long prec;
    Integer val;
    approx_mem.get_best_approximation(361,prec,val);
    assert(prec==38);
    assert(val==Integer(19));

    approx_mem.update_approximation(361,37,15);
    approx_mem.get_best_approximation(361,prec,val);
    assert(prec==38);
    assert(val==Integer(19));

    approx_mem.update_approximation(361,41,1666);
    approx_mem.get_best_approximation(361,prec,val);
    assert(prec==41);
    assert(val==Integer(1666));
  }

  {
    CGAL::CGALi::Best_approximation_cache<String,Integer> approx_mem;

    assert(! approx_mem.approximation_exists("dog"));
    assert(! approx_mem.approximation_exists("cat"));
    assert(! approx_mem.approximation_exists("whatever"));

    approx_mem.update_approximation("duck",12,17);
    approx_mem.update_approximation("elephant",200,12345);

    assert(! approx_mem.approximation_exists("mouse"));
    assert(approx_mem.approximation_exists("duck"));
    long prec;
    Integer val;
    approx_mem.get_best_approximation("elephant",prec,val);
    assert(prec==200);
    assert(val==Integer(12345));
  }

}

int main(int argc,char** argv) {
  
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
