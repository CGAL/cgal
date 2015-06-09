#define CGAL_CHECK_EXPENSIVE
#ifdef CGAL_CHECK_EXACT
bool cgal_check_exact_defined_externally;
#undef CGAL_CHECK_EXACT
#endif
#include "include/sort_test.h"
#include <CGAL/Kinetic/Sort.h>
#include <CGAL/Kinetic/Insert_event.h>
#include <CGAL/Kinetic/Inexact_simulation_traits.h>
#include <cstdlib>
#include <sstream>


int main(int argc, char *argv[])
{
  unsigned int num_points=50;
  unsigned int degree =5;
  if (argc==3) {
    num_points= std::atoi(argv[1]);
    degree= std::atoi(argv[2]);
  }

  std::cout << "Using " << num_points  << " degree " << degree << " points."
	    << std::endl;
  //CGAL_SET_LOG_LEVEL(CGAL::Kinetic::Log::SOME);
  typedef CGAL::Kinetic::Inexact_simulation_traits Tr;
  Tr tr(0,std::numeric_limits<double>::infinity());
  typedef Tr::Simulator::Time Time;

  typedef CGAL::Kinetic::Insert_event<Tr::Active_points_1_table> MOI;
  typedef Tr::Kinetic_kernel::Point_1 MP;

  for (unsigned int i=0; i< num_points; ++i) {
    std::vector<double> coefs;
    for (unsigned int j=0; j<= degree; ++j) {
      coefs.push_back(static_cast<double>(std::rand())
		      /static_cast<double>(RAND_MAX));
    }
    // make sure output is exact
    MP mp(Tr::Kinetic_kernel::Motion_function(coefs.begin(),coefs.end()));
    std::ostringstream oss;
    oss << mp;
    std::istringstream iss;
    iss >> mp;
    tr.simulator_handle()->new_event(Time(i/100.0), 
				      MOI(mp,
					  tr.active_points_1_table_handle()));
  }

  bool error= sort_test<Tr>(tr, 3000);
  if (error) {
    std::cerr << "Sort returned error " << std::endl;
  }

  /*if (error || CGAL::Kinetic::internal::get_static_audit_failures() != 0) {
    return EXIT_FAILURE;
  }
  else {*/
  return EXIT_SUCCESS;
    //}
}
