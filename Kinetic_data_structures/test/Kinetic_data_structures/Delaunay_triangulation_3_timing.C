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
  Simulation_traits::Simulator::Pointer sp= simtr.simulator_pointer();

  KDel kdel(simtr);

  std::ifstream in("data/points_3.n=10,d=3");
  if (!in) {
    std::cerr << "Error opening input file: " << "data/Delaunay_triangulation_3.input" << std::endl;
    return EXIT_FAILURE;
  }
  char buf[1000];
  int nread=0;
  while (true ) {
    in.getline(buf, 1000);
    if (!in) break;
    std::istringstream il(buf);
    Simulation_traits::Kinetic_kernel::Point_3 p;
    il >> p;
    //std::cout << p << std::endl;
    simtr.active_objects_table_pointer()->insert(p);
    ++nread;
  }
  std::cout << "Initializing..." << std::flush;
  kdel.set_has_certificates(true);
  std::cout << "done." << std::endl;

  while (simtr.simulator_pointer()->next_event_time()
	 < simtr.simulator_pointer()->end_time()) {
    sp->set_current_event_number(sp->current_event_number()+1);
    if (sp->current_event_number()%100 == 0) {
      std::cout << "Processed " << sp->current_event_number() << " events." << std::endl;
    }
  }

  /*std::copy(kdel.visitor().begin(), kdel.visitor().end(),
    std::ostream_iterator<std::string>(std::cout, "\n"));*/

  std::cout << "Too late on " << too_late__ << " and filtered " << filtered__ << std::endl;
  std::cout << "Grow " << growth__ << " and shrink " << shrink__ << std::endl; 
  std::cout << "Insert " << queue_insertions__ << " and front " << queue_front_insertions__ << std::endl;
  std::cout << "Sturm created " << sturm_created__ << " and shrink " << sturm_refined__ << std::endl;

  return EXIT_SUCCESS;
};
