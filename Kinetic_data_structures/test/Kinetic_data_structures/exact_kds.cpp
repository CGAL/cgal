#define CGAL_CHECK_EXACTNESS
#define CGAL_CHECK_EXPENSIVE

#include <CGAL/Kinetic/Sort.h>
#include <CGAL/Kinetic/Insert_event.h>
#include <CGAL/Kinetic/Exact_simulation_traits.h>
#ifdef CGAL_USE_CORE
#include <CGAL/Kinetic/CORE_Expr_exact_simulation_traits.h>
#endif
#include <cstdlib>
#include "include/sort_test.h"

int main(int argc, char *argv[]) {
  unsigned int num_points=20;
  unsigned int degree =2;
  if (argc==3) {
    num_points= std::atoi(argv[1]);
    degree= std::atoi(argv[2]);
  }
  std::cout << "Using " << num_points  << " degree " << degree << " points.\n";
  CGAL_SET_LOG_LEVEL(CGAL::Log::NONE);
  bool error=false;
  {
    
    typedef CGAL::Kinetic::Exact_simulation_traits Tr;
    Tr tr(0,1000000);
    typedef Tr::Simulator::Time Time;
    
    typedef CGAL::Kinetic::Insert_event<Tr::Active_points_1_table> MOI;
    typedef Tr::Kinetic_kernel::Point_2 MP;
    typedef  Tr::Kinetic_kernel::Motion_function::NT NT;
    
    for (unsigned int i=0; i< num_points; ++i) {
      std::vector<double> coefs;
      for (unsigned int j=0; j<= degree; ++j) {
	coefs.push_back(static_cast<double>(std::rand())/static_cast<double>(RAND_MAX));
      }
      tr.simulator_handle()->new_event(Time(i/100.0), MOI(MP(Tr::Kinetic_kernel::Motion_function(coefs.begin(), coefs.end()),
							     Tr::Kinetic_kernel::Motion_function(NT(0.0))), 
							  tr.active_points_1_table_handle()));
    }
    error= error || sort_test<Tr>(tr);
  }
#ifdef CGAL_USE_CORE
  {
    
    typedef CGAL::Kinetic::CORE_Expr_exact_simulation_traits Tr;
    Tr tr(0,1000000);
    typedef Tr::Simulator::Time Time;
    
    typedef CGAL::Kinetic::Insert_event<Tr::Active_points_1_table> MOI;
    typedef Tr::Kinetic_kernel::Point_2 MP;
    typedef Tr::Kinetic_kernel::Motion_function::NT NT;
    
    for (unsigned int i=0; i< num_points; ++i) {
      std::vector<double> coefs;
      for (unsigned int j=0; j<= degree; ++j) {
	coefs.push_back(static_cast<double>(std::rand())/static_cast<double>(RAND_MAX));
      }
      tr.simulator_handle()->new_event(Time(i/100.0), MOI(MP(Tr::Kinetic_kernel::Motion_function(coefs.begin(), coefs.end()),
							     Tr::Kinetic_kernel::Motion_function(NT(0.0))), 
							  tr.active_points_1_table_handle()));
    }
    error= error || sort_test<Tr>(tr);
  }
#endif
  
 
  if (error || CGAL::Kinetic::internal::get_static_audit_failures() != 0) {
    return EXIT_FAILURE;
  }
  else {
    return EXIT_SUCCESS;
  }
}
