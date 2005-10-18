#include <CGAL/KDS/basic.h>
#include <limits>
#include <CGAL/Gmpq.h>
#include <CGAL/KDS/Exact_simulation_traits_1.h>


template <class NT>
struct Check{
  Check(NT v){
    NT zero(0);
    assert(v <= 0 || v > 0); // break symmetry due to compiler warning
    if (std::numeric_limits<NT>::has_infinity){
      NT inf= std::numeric_limits<NT>::infinity();
      NT inf2=inf;
      NT ninf= -inf;
      assert(inf > v);
      assert(ninf < v);
      assert(inf > zero);
      assert(ninf < zero);
      assert(inf > ninf);
      assert(ninf < inf);
      assert(inf2==inf);
    }
  }
};

int main(int, char *[]){
  Check<double> d(0);
  Check<CGAL::Gmpq> g(1);

  typedef CGAL::KDS::Exact_simulation_traits_1 Tr;
  Tr tr;
  Tr::Simulator::Function_kernel::Function fn= tr.simulator_pointer()->function_kernel_object().construct_function_object()(1,0,-2);
  Check<Tr::Simulator::Time> r(tr.simulator_pointer()->root_stack_object(fn).top());


  return EXIT_SUCCESS;
}
