#define CGAL_CHECK_EXPENSIVE
#define CGAL_CHECK_EXACTNESS

#include <CGAL/IO/Qt_debug_viewer_2.h>
#include <CGAL/Arrangement_of_spheres_3/Slice.h>
#include <CGAL/Arrangement_of_spheres_traits_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Gmpq.h>
#include <iostream>
#include <sstream>


int main(int argc, char *argv[]){

  //typedef CGAL::Exact_predicates_exact_constructions_kernel K;
  typedef Arrangement_of_spheres_traits_3::Geometric_traits K;
  typedef CGAL::Simple_cartesian<double> DK;
  std::vector<K::Sphere_3> spheres;
  std::vector<K::Sphere_3> shifted_spheres;
  while (true){
    char buf[1000];
    std::cin.getline(buf, 1000);
    if (!std::cin) break;
    std::istringstream iss(buf);
    DK::Sphere_3 s;
    iss >> s;
    if (!iss) {
      std::cerr << "Can't parse line " << buf << std::endl;
    } else {
      spheres.push_back(K::Sphere_3(K::Point_3(s.center().x(), s.center().y(),
					       s.center().z()), s.squared_radius()));
      shifted_spheres.push_back(K::Sphere_3(K::Point_3(s.center().x()+.1, s.center().y()+.1,
						       s.center().z()), s.squared_radius()));
    }
    //std::cout << spheres.back() << std::endl;
  }


  std::cout << "Read " << spheres.size() << " spheres." << std::endl;
  Slice::T tr(spheres.begin(), spheres.end());
  std::cout << "Bounding box is from " << tr.bbox_3().xmin() << " to " 
	    << tr.bbox_3().xmax()
	    << std::endl;
 
  Geometric_traits::FT z= atof(argv[1]);
 
  //
 
  Slice slice(tr);
  slice.initialize_at(z);

  QApplication app(argc, argv);
  Qt_examiner_viewer_2 *qtd= new Qt_examiner_viewer_2(10);
 
  
   
  slice.draw_rz(qtd, z);
  app.setMainWidget( qtd);
  qtd->show_everything();
  qtd->show();
  return app.exec();
}
