#include <CGAL/KDS/basic.h>

#include <CGAL/KDS/Delaunay_triangulation_3.h>
#include <CGAL/KDS/Exact_simulation_traits_3.h>
#include <CGAL/KDS/Delaunay_triangulation_event_log_visitor_3.h>
#include <algorithm>
#include <iterator>

int main(int, char *[]){

  typedef CGAL::KDS::Exact_simulation_traits_3 Simulation_traits;
  typedef CGAL::KDS::Delaunay_triangulation_event_log_visitor_3 Visitor;
  typedef CGAL::KDS::Delaunay_triangulation_3<Simulation_traits, Visitor> KDel;


  Simulation_traits simtr;
  Simulation_traits::Simulator::Pointer sp= simtr.simulator_pointer();

  KDel kdel(simtr);

  std::ifstream in("data/Delaunay_triangulation_3.input");
  if (!in) {
    std::cerr << "Error opening input file: " << "data/Delaunay_triangulation_3.input" << std::endl;
    return EXIT_FAILURE;
  }
  char buf[1000];
  int nread=0;
  while (true ){
    in.getline(buf, 1000);
    if (!in) break;
    std::istringstream il(buf);
    Simulation_traits::Kinetic_kernel::Point_3 p;
    il >> p;
    //std::cout << p << std::endl;
    simtr.moving_point_table_pointer()->insert(p);
    ++nread;
  }
  kdel.set_has_certificates(true);

  while (simtr.simulator_pointer()->next_event_time() 
	 < simtr.simulator_pointer()->end_time()){
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
  for (CGAL::KDS::Delaunay_triangulation_event_log_visitor_3::iterator it = kdel.visitor().begin();
       it != kdel.visitor().end(); ++it){
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
  

  if (error_count != 0){
    std::cerr << "ERROR " << error_count << " errors in " << kdel.visitor().size() << " events.\n";
  }
  return EXIT_SUCCESS;
};
