#include <CGAL/KDS/Exact_simulation_traits_1.h>
#include <CGAL/KDS/Inexact_simulation_traits_1.h>

template <class Traits, class Fn, class Rt>
void check_one(const Traits &tr, const Fn &fn, const Rt &lb, const Rt &rt) {
  typename Traits::Root_stack rs= tr.root_stack_object(fn, lb, std::numeric_traits<Rt>::infinity());
  
  if (rt != std::numeric_traits<Rt>::infinity()) {
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
    CGAL::KDS::Exact_simulation_traits_1 tr;
    typedef CGAL::KDS::Exact_simulation_traits_1::Root Root;
    Root inf= std::numeric_limits<Root>::infinity();
    Traits::Construct_function cf= tr.construct_function_object();
    check_one(tr, cf(1,-2,1), 0, inf);
  }

  {
    typedef CGAL::KDS::Inexact_simulation_traits_1 Traits;
  }

  return EXIT_SUCCESS;
}
