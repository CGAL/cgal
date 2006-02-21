//#define CGAL_CHECK_EXPENSIVE
//#define CGAL_CHECK_EXACTNESS

#include <CGAL/Kinetic/basic.h>

#include <CGAL/Kinetic/Delaunay_triangulation_3.h>
#include <CGAL/Kinetic/Exact_simulation_traits_3.h>
#include <algorithm>
#include <iterator>

int main(int, char *[])
{

  typedef CGAL::Kinetic::Exact_simulation_traits_3 Simulation_traits;
  typedef CGAL::Kinetic::Delaunay_triangulation_3<Simulation_traits> KDel;

  Simulation_traits simtr;
  Simulation_traits::Simulator::Pointer sp= simtr.simulator_pointer();

  KDel kdel(simtr);

  std::ifstream in("data/points_3.n=100,d=1");
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
  kdel.set_has_certificates(true);

  while (simtr.simulator_pointer()->next_event_time()
	 < simtr.simulator_pointer()->end_time()) {
    sp->set_current_event_number(sp->current_event_number()+1);
  }

  /*std::copy(kdel.visitor().begin(), kdel.visitor().end(),
    std::ostream_iterator<std::string>(std::cout, "\n"));*/

  std::cout << "Too late on " << too_late__ << " and filtered " << filtered__ << std::endl;

  return EXIT_SUCCESS;
};
