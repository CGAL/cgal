#define CGAL_CHECK_EXPENSIVE
#ifdef CGAL_CHECK_EXACT
bool cgal_check_exact_defined_externally;
#undef CGAL_CHECK_EXACT
#endif

#include <CGAL/KDS/Sort.h>
#include <CGAL/KDS/Insert_event.h>
#include <CGAL/KDS/Inexact_simulation_traits_1.h>
#include <cstdlib>
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
    std::cout << "Using " << num_points  << " degree " << degree << " points.\n";
    //CGAL_KDS_SET_LOG_LEVEL(CGAL::KDS::LOG_SOME);
    typedef CGAL::KDS::Inexact_simulation_traits_1 Tr;
    Tr tr;
    typedef Tr::Simulator::Time Time;

    typedef CGAL::KDS::Insert_event<Tr::Active_objects_table> MOI;
    typedef Tr::Kinetic_kernel::Point_1 MP;

    for (unsigned int i=0; i< num_points; ++i) {
      std::vector<double> coefs;
      for (unsigned int j=0; j<= degree; ++j) {
	coefs.push_back(static_cast<double>(std::rand())
			/static_cast<double>(RAND_MAX));
      }
      tr.simulator_pointer()->new_event(Time(i/100.0), 
					MOI(MP(Tr::Function_kernel::Function(coefs.begin(),coefs.end())),
					    tr.active_objects_table_pointer()));
    }

    error= sort_test<Tr>(tr);
  } while (false);


  if (error) {
    return EXIT_FAILURE;
  }
  else {
    return EXIT_SUCCESS;
  }
};
