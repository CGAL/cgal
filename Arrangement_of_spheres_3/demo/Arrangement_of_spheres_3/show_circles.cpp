#define CGAL_CHECK_EXPENSIVE
#define CGAL_CHECK_EXACTNESS

#include <CGAL/IO/Qt_debug_viewer_2.h>
#include <CGAL/Arrangement_of_spheres_3/Slice_arrangement.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arrangement_of_spheres_3/coordinates.h>
#include <iostream>
#include <sstream>

typedef CGAL::Simple_cartesian<double> K;


double squared_depth(K::Point_3 p, double z){
  return CGAL::square(p[sweep_coordinate().index()]-z);
}

bool has_overlap(const K::Sphere_3 &s,  double z) {
  double dz2= squared_depth(s.center(), z);
  if (dz2 <= s.squared_radius()) return true;
  else return false;
}

K::Circle_2 intersect(K::Sphere_3 s, double z){
  double r2= s.squared_radius() - squared_depth(s.center(), z);
  CGAL_assertion(r2>=0);
  K::Circle_2  c(K::Point_2(s.center()[plane_coordinate(0).index()],
			    s.center()[plane_coordinate(1).index()]), r2);
  return c;
}


int main(int argc, char *argv[]){
 
  
  CGAL::Bbox_3 box(std::numeric_limits<double>::max(),
		   std::numeric_limits<double>::max(),
		   std::numeric_limits<double>::max(),
		   -std::numeric_limits<double>::max(),
		   -std::numeric_limits<double>::max(),
		   -std::numeric_limits<double>::max());
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
      box= box+ spheres.back().bbox();
    }
    //std::cout << spheres.back() << std::endl;
  }


  std::cout << "Read " << spheres.size() << " spheres." << std::endl;
  std::cout << "Bounding box is from " << box.xmin() << " to " << box.xmax()
	    << std::endl;
 
  Geometric_traits::FT z= atof(argv[1]);


  QApplication app(argc, argv);
  Qt_examiner_viewer_2 *qtd= new Qt_examiner_viewer_2(10);
 
  for (unsigned int i=0; i< spheres.size(); ++i){
    if (has_overlap(spheres[i], z)) {
      K::Circle_2 c= intersect(spheres[i], z);
      *qtd << c;
    }
  }
   
 
  app.setMainWidget( qtd);
  qtd->show_everything();
  qtd->show();
  return app.exec();
}
