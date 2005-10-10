#define NDEBUG

#include <CGAL/KDS/Exact_simulation_traits_1.h>
#include <CGAL/KDS/Inexact_simulation_traits_1.h>

template <class Traits, class Fn, class Rt>
void check_one(const Traits &tr, const Fn &fn, const Rt &lb, const Rt &rt) {
  typename Traits::Simulator::Root_stack rs(fn, lb, std::numeric_limits<Rt>::infinity(), 
					    tr.simulator_pointer()->function_kernel_object());
  
  if (rt != std::numeric_limits<Rt>::infinity()) {
    if (rs.top() != rt) {
      std::cerr << "ERROR For function " << fn << " expected " << rt << " got " << rs.top() << std::endl;
    }
  } else {
    if (!rs.empty()){
      std::cerr << "ERROR For function " << fn << " expected " << rt << " got " << rs.top() << std::endl;
    }
  }
}


//check the KDS solvers

int main(int, char *[]){


  {
    typedef CGAL::KDS::Exact_simulation_traits_1 Traits;
    Traits tr;
    typedef Traits::Simulator::Root_stack::Root Root;
    Root inf= std::numeric_limits<Root>::infinity();
    Traits::Simulator::Function_kernel::Construct_function cf= tr.function_kernel_object().construct_function_object();
    check_one(tr, cf(1,-2,1), Root(0), inf);
  }

  {
    typedef CGAL::KDS::Inexact_simulation_traits_1 Traits;
  }

  return EXIT_SUCCESS;
}
