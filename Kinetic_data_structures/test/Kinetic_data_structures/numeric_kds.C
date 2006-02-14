#define CGAL_CHECK_EXPENSIVE
#ifdef CGAL_CHECK_EXACT
bool cgal_check_exact_defined_externally;
#undef CGAL_CHECK_EXACT
#endif

#include <CGAL/Kinetic/Sort.h>
#include <CGAL/Kinetic/Insert_event.h>
#include <CGAL/Kinetic/Inexact_simulation_traits_1.h>
#include <cstdlib>
#include <sstream>
#include "include/sort_test.h"

int main(int argc, char *argv[])
{
  unsigned int num_points=50;
  unsigned int degree =5;
  /*if (argc==3) {
    num_points= std::atoi(argv[1]);
    degree= std::atoi(argv[2]);
  }*/

  bool error;

  do {
    std::cout << "Using " << num_points  << " degree " << degree << " points."
	<< std::endl;
    //CGAL_KINETIC_SET_LOG_LEVEL(CGAL::Kinetic::LOG_SOME);
    typedef CGAL::Kinetic::Inexact_simulation_traits_1 Tr;
    Tr tr;
    typedef Tr::Simulator::Time Time;

    typedef CGAL::Kinetic::Insert_event<Tr::Active_objects_table> MOI;
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
      tr.simulator_pointer()->new_event(Time(i/100.0), 
					MOI(mp,
					    tr.active_objects_table_pointer()));
    }

    error= sort_test<Tr>(tr, 3000);
  } while (false);


  if (error) {
    return EXIT_FAILURE;
  }
  else {
    return EXIT_SUCCESS;
  }
};
