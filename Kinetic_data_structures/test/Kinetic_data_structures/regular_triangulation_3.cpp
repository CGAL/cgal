#define CGAL_CHECK_EXACTNESS
#define CGAL_CHECK_EXPENSIVE

#include <CGAL/Kinetic/basic.h>

#include <CGAL/Kinetic/Regular_triangulation_3.h>
#include <CGAL/Kinetic/Regular_triangulation_CORE_exact_simulation_traits.h>
#include <CGAL/Kinetic/Regular_triangulation_inexact_simulation_traits.h>
#include <CGAL/Kinetic/Regular_triangulation_event_log_visitor_3.h>
#include <algorithm>
#include <iterator>

int main(int, char *[])
{

 

  if (0) {
    typedef CGAL::Kinetic::Regular_triangulation_inexact_simulation_traits Simulation_traits;
    typedef CGAL::Kinetic::Regular_triangulation_event_log_visitor_3 Visitor;
    typedef CGAL::Kinetic::Regular_triangulation_3<Simulation_traits, Visitor> KDel;
    Simulation_traits simtr(0,100000);
    Simulation_traits::Simulator::Handle sp= simtr.simulator_handle();
      
    KDel kdel(simtr);
    //CGAL_KINETIC_SET_LOG_LEVEL(CGAL::Kinetic::LOG_LOTS);
    std::ifstream in("data/regular_triangulation_3.input");
    if (!in) {
      std::cerr << "Error opening input file: " << "data/regular_triangulation_3.input" << std::endl;
      return EXIT_FAILURE;
    }
    char buf[1000];
    int nread=0;
    while (true ) {
      in.getline(buf, 1000);
      if (!in) break;
      std::istringstream il(buf);
      Simulation_traits::Kinetic_kernel::Weighted_point_3 p;
      il >> p;
      //std::cout << p << std::endl;
      simtr.active_points_3_table_handle()->insert(p); // here 
      ++nread;
    }
    kdel.set_has_certificates(true);
      
    while (sp->next_event_time()
	   < sp->end_time()) {
      sp->set_current_event_number(sp->current_event_number()+1);
    }
  } else {
    typedef CGAL::Kinetic::Regular_triangulation_CORE_exact_simulation_traits Simulation_traits;
    typedef CGAL::Kinetic::Regular_triangulation_event_log_visitor_3 Visitor;
    typedef CGAL::Kinetic::Regular_triangulation_3<Simulation_traits, Visitor> KDel;
    Simulation_traits simtr(0,100000);
    Simulation_traits::Simulator::Handle sp= simtr.simulator_handle();

    KDel kdel(simtr);
    //CGAL_KINETIC_SET_LOG_LEVEL(CGAL::Kinetic::LOG_LOTS);
    std::ifstream in("data/regular_triangulation_3.input");
    if (!in) {
      std::cerr << "Error opening input file: " << "data/regular_triangulation_3.input" << std::endl;
      return EXIT_FAILURE;
    }
    char buf[1000];
    int nread=0;
    while (true ) {
      in.getline(buf, 1000);
      if (!in) break;
      std::istringstream il(buf);
      Simulation_traits::Kinetic_kernel::Weighted_point_3 p;
      il >> p;
      //std::cout << p << std::endl;
      simtr.active_points_3_table_handle()->insert(p); // here 
      ++nread;
    }
    kdel.set_has_certificates(true);

    while (sp->next_event_time()
	   < sp->end_time()) {
      sp->set_current_event_number(sp->current_event_number()+1);
    }

    /*std::copy(kdel.visitor().begin(), kdel.visitor().end(),
      std::ostream_iterator<std::string>(std::cout, "\n"));*/

    std::ifstream out("data/regular_triangulation_3.output");
    if (!out) {
      std::cerr << "Error opening input file: " << "data/regular_triangulation_3.output" << std::endl;
      return EXIT_FAILURE;
    }

    int error_count=0;
    for (CGAL::Kinetic::Delaunay_triangulation_event_log_visitor_3::Event_iterator it = kdel.visitor().events_begin();
	 it != kdel.visitor().events_end(); ++it) {
      char buf[1000];
      out.getline(buf, 1000);
      if (*it != buf) {
	std::cerr << "ERROR Got event: " << *it << std::endl;
	std::cerr << "      Expected event: " << buf << std::endl;
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
    }
    CGAL::Kinetic::internal::write_debug_counters(std::cout);
    if (CGAL::Kinetic::internal::audit_failures__ != 0) return EXIT_FAILURE;
    else return EXIT_SUCCESS;
  }
};
