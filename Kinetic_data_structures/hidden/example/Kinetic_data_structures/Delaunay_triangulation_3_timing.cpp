//#define CGAL_CHECK_EXPENSIVE
//#define CGAL_CHECK_EXACTNESS

#include <CGAL/Kinetic/basic.h>

#include <CGAL/Kinetic/Delaunay_triangulation_3.h>
#include <CGAL/Kinetic/Exact_simulation_traits_3.h>
#include <CGAL/Kinetic/Inexact_simulation_traits_3.h>
#include <algorithm>
#include <iterator>

int main(int, char *[])
{
  typedef CGAL::Kinetic::Exact_simulation_traits_3 Simulation_traits;
  //typedef CGAL::Kinetic::Inexact_simulation_traits_3 Simulation_traits;
  typedef CGAL::Kinetic::Delaunay_triangulation_3<Simulation_traits> KDel;

  Simulation_traits simtr;
  Simulation_traits::Simulator::Handle sp= simtr.simulator_handle();

  KDel kdel(simtr);

  std::string inputf("data/points_3.n=100,d=1");

  std::ifstream in(inputf.c_str());
  if (!in) {
    std::cerr << "Error opening input file: " << inputf << std::endl;
    return EXIT_FAILURE;
  }

  in >> *simtr.active_points_3_table_handle();
  //std::cout <<  *simtr.active_objects_table_pointer();
  
  std::cout << "Read " <<  simtr.active_points_3_table_handle()->size() << std::endl;
  std::cout << "Initializing..." << std::flush;
  kdel.set_has_certificates(true);
  std::cout << "done." << std::endl;

  while (simtr.simulator_handle()->next_event_time()
	 < simtr.simulator_handle()->end_time()) {
    sp->set_current_event_number(sp->current_event_number()+1);
    if (sp->current_event_number()%100 == 0) {
      std::cout << "Processed " << sp->current_event_number() << " events." << std::endl;
    }
  }

  std::cout << "Proccessed " << sp->current_event_number() << " events" << std::endl;
  /*std::copy(kdel.visitor().begin(), kdel.visitor().end(),
    std::ostream_iterator<std::string>(std::cout, "\n"));*/

  //std::cout << "Too late on " << too_late__ << " and filtered " << filtered__ << std::endl;
  //  std::cout << "Grow " << growth__ << " and shrink " << shrink__ << std::endl; 
  //std::cout << "Insert " << queue_insertions__ << " and front " << queue_front_insertions__ << std::endl;
  //std::cout << "Sturm created " << sturm_created__ << " and shrink " << sturm_refined__ << std::endl;

  
  if (CGAL::Kinetic::internal::fail__) return EXIT_FAILURE;
  else return EXIT_SUCCESS;
};
