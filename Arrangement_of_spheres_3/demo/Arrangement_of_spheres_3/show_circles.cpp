//#define CGAL_CHECK_EXPENSIVE
//#define CGAL_CHECK_EXACTNESS
#include <CGAL/Arrangement_of_spheres_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Qt_multithreaded_examiner_viewer_2.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_qt_viewer.h>
#include <CGAL/Arrangement_of_spheres_3/read_spheres.h>
#include <vector>

typedef CGAL::Simple_cartesian<double> K;


double squared_depth(K::Point_3 p, double x){
  return CGAL::square(p[0]-x);
}

bool has_overlap(const K::Sphere_3 &s,  double x) {
  double dx2= squared_depth(s.center(), x);
  if (dx2 <= s.squared_radius()) return true;
  else return false;
}

K::Circle_2 intersect(K::Sphere_3 s, double x){
  double r2= s.squared_radius() - squared_depth(s.center(), x);
  CGAL_assertion(r2>=0);
  K::Circle_2  c(K::Point_2(s.center()[1],
			    s.center()[2]), r2);
  return c;
}


struct Show_circles {
  Show_circles(const std::vector<K::Sphere_3> &ss, double z):
    spheres_(ss), z_(z){}

  void operator()(CGAL::Qt_examiner_viewer_2 *qtv) {
  
    for (unsigned int i=0; i< spheres_.size(); ++i) {
      if (has_overlap(spheres_[i], z_)) {
	*qtv << intersect(spheres_[i], z_);
	std::cout << intersect(spheres_[i], z_) << std::endl;
      }
    }

    qtv->show_everything();
  }

  std::vector< K::Sphere_3> spheres_;
  double z_;
};

 

int main(int argc, char *argv[]) {
#ifdef CGAL_AOS3_USE_TEMPLATES
  //typedef CGAL::Simple_cartesian<double> K;
  typedef CGAL::Arrangement_of_spheres_traits_3<K> Traits;
  typedef CGAL::Arrangement_of_spheres_3<Traits> Arrangement;
#else 
  typedef CGAL::Arrangement_of_spheres_3 Arrangement;
#endif

  std::vector<K::Sphere_3> spheres;

  CGAL_AOS3_INTERNAL_NS::read_spheres<K, true>(std::cin, spheres);
  
  double z=std::atof(argv[1]);

  typedef Show_circles CS;
  CS cs(spheres, z);
  
  CGAL::Qt_multithreaded_examiner_viewer_2<CS> qtv(cs, argc, argv);
  

  return qtv();
}
