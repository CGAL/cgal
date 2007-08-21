#include <CGAL/Arrangement_of_spheres_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Qt_multithreaded_examiner_viewer_2.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_qt_viewer.h>
#include <CGAL/Arrangement_of_spheres_3/read_spheres.h>
#include <vector>

typedef CGAL::Simple_cartesian<double> DK;

#ifdef CGAL_AOS3_USE_TEMPLATES
  typedef CGAL::Simple_cartesian<double> K;
  typedef CGAL::Arrangement_of_spheres_traits_3<K> Traits;
  typedef CGAL::Arrangement_of_spheres_3<Traits> Arrangement;
#else 
  typedef CGAL::Arrangement_of_spheres_3 Arrangement;
#endif

template <class A>
struct Compute_slice {
  Compute_slice(const std::vector<typename Arrangement::Traits::Geom_traits::Sphere_3> &ss, double z):
    spheres_(ss), z_(z){}

  void operator()(CGAL::Qt_examiner_viewer_2 *qtv) {
    typedef typename A::Traits::Geom_traits K;
    typedef typename A::Traits Traits_t;
    typedef typename K::Point_3 P;
    typedef typename K::Sphere_3 S;
   
    std::vector<S> ns;
    for (unsigned int i=0; i< spheres_.size(); ++i) {
      ns.push_back(S(P(spheres_[i].center().x(),
		       spheres_[i].center().y(),
		       spheres_[i].center().z()),
		     spheres_[i].squared_radius()));
      std::cout << ns.back() << std::endl;
    }

    typename A::Traits tr(ns.begin(), ns.end());

    typename A::Cross_section cs;
   
    cs.visitor().set_start_time(end);
    cs.visitor().set_traits(tr_);
    cs.visitor().set_cross_section(cs);
    cs.set_number_of_spheres(tr_.number_of_sphere_3s());
    std::cout << "Initialize at " << end << std::endl;
    CGAL_AOS3_INTERNAL_NS::Cross_section_initializer CGAL_AOS3_TARG csi(cs, tr_);
    csi(end, qtv);
  }

  std::vector<typename Arrangement::Traits::Geom_traits::Sphere_3> spheres_;
  double z_;
};

 

int main(int argc, char *argv[]) {


  std::vector<Arrangement::Traits::Geom_traits::Sphere_3> spheres;

  CGAL_AOS3_INTERNAL_NS::read_spheres<Arrangement::Traits::Geom_traits, true>(std::cin, spheres);
  
  double z=std::atof(argv[1]);

  typedef Compute_slice<Arrangement> CS;
  CS cs(spheres, z);
  
  CGAL::Qt_multithreaded_examiner_viewer_2<CS> qtv(cs, argc, argv);
  

  return qtv();
}
