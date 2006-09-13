#define CGAL_CHECK_EXPENSIVE
#define CGAL_CHECK_EXACTNESS

#include <CGAL/IO/Qt_debug_viewer_2.h>
#include <CGAL/Arrangement_of_spheres_3/Slice_arrangement.h>
#include <CGAL/Arrangement_of_spheres_traits_3.h>
#include <CGAL/Arrangement_of_spheres_3/Slice.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Gmpq.h>
#include <iostream>
#include <sstream>
#include <fstream>





int main(int argc, char *argv[]){
  
  typedef Arrangement_of_spheres_traits_3::Geometric_traits K;
  typedef CGAL::Simple_cartesian<double> DK;
  std::vector<K::Sphere_3> spheres;
  std::ifstream in(argv[1]);
  while (true){
    char buf[1000];
    in.getline(buf, 1000);
    if (!in) break;
    std::istringstream iss(buf);
    DK::Sphere_3 s;
    iss >> s;
    if (!iss) {
      std::cerr << "Can't parse line " << buf << std::endl;
    } else {
      spheres.push_back(K::Sphere_3(K::Point_3(s.center().x(), 
					       s.center().y(),
					       s.center().z()), 
				    s.squared_radius()));
    }
    //std::cout << spheres.back() << std::endl;
  }



  std::cout << "Read " << spheres.size() << " spheres." << std::endl;
  Arrangement_of_spheres_traits_3 tr(spheres.begin(), spheres.end());

  std::cout << "Bounding box is from " << tr.bbox_3().xmin() << " to " << tr.bbox_3().xmax()
	    << std::endl;

  //Simulator::Handle s= new Simulator(Simulator::Time(tr.bbox_3().xmin()-1), Simulator::Time(tr.bbox_3().xmax()+1));
 
  

  
  Slice slice(tr);
  Qt_gui::Handle qtg= new Qt_gui(argc, argv, slice.simulator_handle());
  slice.set_gui(qtg);
  
  return qtg->begin_event_loop();
}
