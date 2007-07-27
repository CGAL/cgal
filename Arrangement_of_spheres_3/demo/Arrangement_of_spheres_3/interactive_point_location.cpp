#define CGAL_CHECK_EXPENSIVE
#define CGAL_CHECK_EXACTNESS

#include <CGAL/Arrangement_of_spheres_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Qt_multithreaded_examiner_viewer_2.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_qt_viewer.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_qt_viewer_markup.h>
#include <CGAL/Arrangement_of_spheres_3/read_spheres.h>
#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section.h>
#include <vector>


template <class A>
struct Do_work {
  Do_work(const std::vector<typename A::Traits::Geom_traits::Sphere_3> &ss, double z):
    spheres_(ss), z_(z){}

  void operator()(CGAL::Qt_examiner_viewer_2 *q){
    typedef typename A::Traits::Geom_traits K;
    typedef typename A::Traits Traits_t;
    typedef typename K::Point_3 P;
    typedef typename K::Sphere_3 S;

    CGAL_AOS3_INTERNAL_NS::Unprojector<K> up(z_);
   
    std::vector<S> ns;
    for (unsigned int i=0; i< spheres_.size(); ++i) {
      ns.push_back(S(P(spheres_[i].center().x(),
		       spheres_[i].center().y(),
		       spheres_[i].center().z()),
		     spheres_[i].squared_radius()));
      std::cout << ns.back() << std::endl;
    }

    typename A::Traits tr(ns.begin(), ns.end());
    A arr(tr);
    typename A::Cross_section cs(tr.number_of_sphere_3s());
    arr.initialize_at(z_,cs);

    typedef typename CGAL_AOS3_INTERNAL_NS::Cross_section_qt_viewer CGAL_AOS3_TARG CSV;
    CSV csv(tr, cs);

  
    *q << CGAL::Layer(0);
    csv(z_, q);
    
    q->show_everything();
    //q->redraw();
 

    typedef typename CGAL_AOS3_INTERNAL_NS::Irrational_cross_section CGAL_AOS3_TARG ICS;
    ICS ics(tr, cs, z_);

    typedef typename CGAL_AOS3_INTERNAL_NS::Cross_section_qt_viewer_markup CGAL_AOS3_TARG MCSV;
    MCSV mcsv(tr, cs, CGAL::Layer(1));

    while (true) {
      std::cout << "Enter coordinates: " << std::flush;
      char buf[1000];
      std::cin.getline(buf,1000);
      if (buf[0]== '\0') {
	std::cout << "bye." << std::endl;
	break;
      }
      mcsv.clear();

      double x,y;
      std::istringstream iss(buf);
      iss >> x >> y;
      if (!iss) {
	std::cerr << "Can't parse line." << std::endl;
      } else {
	*q << CGAL::Layer(2);
	*q << CGAL::RED;
	*q << typename K::Point_2(x,y);
	//q->redraw();
	typename A::Traits::Sphere_point_3 sp(up(typename K::Point_2(x,y)), 
					      typename K::Line_3(up(typename K::Point_2(x,y)),
								 CGAL_AOS3_INTERNAL_NS::Sweep_coordinate::basis_3<typename K::Vector_3>()));
	try {
	  typename A::Cross_section::Face_const_handle f= ics.locate_point(sp);
	  //slice.new_marked_face(f);
	  std::cout << "Found face ";
	  cs.write(f, std::cout) << std::endl;
	  mcsv.new_face(f);

	} catch (ICS::On_edge_exception e) {
	  std::cout << "On edge!" <<std::endl;
	  mcsv.new_edge(e.halfedge_handle());
	} catch (ICS::On_vertex_exception v) {
	  std::cout << "On vertex!" <<std::endl;
	  mcsv.new_vertex(v.vertex_handle());
	}
	mcsv(z_, q);
      }
    }
  }
  

  std::vector<typename A::Traits::Geom_traits::Sphere_3> spheres_;
  double z_;

};


int main(int argc, char *argv[]){
#ifdef CGAL_AOS3_USE_TEMPLATES
  typedef CGAL::Simple_cartesian<double> K;
  typedef CGAL::Arrangement_of_spheres_traits_3<K> Traits;
  typedef CGAL::Arrangement_of_spheres_3<Traits> Arrangement;
#else 
  typedef CGAL::Arrangement_of_spheres_3 Arrangement;
#endif

  std::vector<Arrangement::Traits::Geom_traits::Sphere_3> spheres;

  std::ifstream in(argv[2]);

  CGAL_AOS3_INTERNAL_NS::read_spheres<Arrangement::Traits::Geom_traits, true>(in, spheres);
  
  double z=std::atof(argv[1]);

  typedef Do_work<Arrangement> CS;
  CS cs(spheres, z);
  
  CGAL::Qt_multithreaded_examiner_viewer_2<CS> qtv(cs, argc, argv);
  

  return qtv();
}
