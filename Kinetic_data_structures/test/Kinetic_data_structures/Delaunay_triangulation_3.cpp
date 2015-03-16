#define CGAL_CHECK_EXACTNESS
#define CGAL_CHECK_EXPENSIVE

#include <CGAL/Kinetic/basic.h>

#include <CGAL/Kinetic/Delaunay_triangulation_3.h>
#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/Delaunay_triangulation_event_log_visitor_3.h>
#include <algorithm>
#include <iterator>

int main(int argc, char *argv[])
{

  typedef CGAL::Kinetic::Exact_simulation_traits Simulation_traits;
  typedef CGAL::Kinetic::Delaunay_triangulation_event_log_visitor_3 Visitor;
  typedef CGAL::Kinetic::Delaunay_triangulation_3<Simulation_traits, Visitor> KDel;

  Simulation_traits simtr(0,100000);
  Simulation_traits::Simulator::Handle sp= simtr.simulator_handle();

  KDel kdel(simtr);
  if (argc==1) {
    CGAL_SET_LOG_LEVEL(CGAL::Log::NONE);
    std::ifstream in("data/Delaunay_triangulation_3.input");
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
      simtr.active_points_3_table_handle()->insert(p);
      ++nread;
    }
    kdel.set_has_certificates(true);
    kdel.audit();

  } else {
    CGAL_SET_LOG_LEVEL(CGAL::Log::LOTS);
    std::ifstream in(argv[1]);
    if (!in) {
      std::cerr << "Error opening input file: " << argv[1] << std::endl;
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
      simtr.active_points_3_table_handle()->insert(p);
      ++nread;
    }
    kdel.set_has_certificates(true);
    kdel.audit();
  }

  while (simtr.simulator_handle()->next_event_time()
	 < simtr.simulator_handle()->end_time()) {
    sp->set_current_event_number(sp->current_event_number()+1);
  }

  /*std::copy(kdel.visitor().begin(), kdel.visitor().end(),
    std::ostream_iterator<std::string>(std::cout, "\n"));*/

  std::ifstream out("data/Delaunay_triangulation_3.output");
  if (!out) {
    std::cerr << "Error opening input file: " << "data/Delaunay_triangulation_3.output" << std::endl;
    return EXIT_FAILURE;
  }

  int error_count=0;
  for (CGAL::Kinetic::Delaunay_triangulation_event_log_visitor_3::Event_iterator it = kdel.visitor().events_begin();
       it != kdel.visitor().events_end(); ++it) {
    char buf[1000];
    out.getline(buf, 1000);
    if (*it != buf) {
      std::cerr << "ERROR Got event: " << *it << std::endl;
      std::cerr << "ERROR Expected event: " << buf << std::endl;
      ++error_count;
    }
  }

  while (out) {
    char buf[1000];
    out.getline(buf, 1000);
    if (out) {
      std::cerr << "ERROR Missing event: " << buf << std::endl;
      ++error_count;
    }
  }

  if (error_count != 0) {
    std::cerr << "ERROR " << error_count << " errors in " << kdel.visitor().size() << " events.\n";
  } else {
    std::cout << "No errors for " << kdel.visitor().size() << " events" << std::endl;
  }

  //std::cout << "Too late on " << too_late__ << " and filtered " << filtered__ << std::endl;
  //std::cout << "Grow " << growth__ << " and shrink " << shrink__ << std::endl;
  //std::cout << "Insert " << queue_insertions__ << " and front " << queue_front_insertions__ << std::endl;
  //std::cout << "Sturm created " << sturm_created__ << " and shrink " << sturm_refined__ << std::endl;
  
  if (CGAL::Kinetic::internal::get_static_audit_failures() != 0 ) return EXIT_FAILURE;
  else return EXIT_SUCCESS;
}
