#define CGAL_CHECK_EXPENSIVE
#define CGAL_CHECK_EXACTNESS

#include <CGAL/Arrangement_of_spheres_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Qt_multithreaded_examiner_viewer_2.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_qt_viewer.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_qt_viewer_markup.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_qt_viewer_events.h>
#include <CGAL/Arrangement_of_spheres_3/read_spheres.h>
#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section_location.h>
#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section_insertion.h>
#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section_removal.h>
#include <vector>


void my_failure_handler(
			const char */*type*/,
			const char *expr,
			const char* file,
			int line,
			const char* msg)
{
  std::cerr << "ASSERTION FAILURE " << expr << " in file " << file 
	    << ": " << line << std::endl << msg << std::endl;
  //int i;

  while (true){};
}


template <class T, class CS>
void dump(const T&tr, const CS &cs) {
  for (typename T::Sphere_3_key_const_iterator it= tr.sphere_3_keys_begin();
       it != tr.sphere_3_keys_end(); ++it) {
    if (cs.a_halfedge(*it) != typename CS::Halfedge_handle()) {
      std::cout << tr.sphere_3(*it) << std::endl;
    }
  }
}

template <class A>
struct Do_work {
  Do_work(const std::vector<typename A::Traits::Geom_traits::Sphere_3> &ss, double z):
    spheres_(ss), z_(z){}

  void operator()(CGAL::Qt_examiner_viewer_2 *q){
    typedef typename A::Traits::Geom_traits K;
    typedef typename A::Traits Traits_t;
    typedef typename Traits_t::Event_key Event_key;
    typedef typename K::Point_3 P;
    typedef typename K::Sphere_3 S;

    double zcur= CGAL::to_double(z_);

    CGAL_AOS3_INTERNAL_NS::Unprojector<K> up(z_);
    //CGAL_AOS3_INTERNAL_NS::Projector<K> down(z_);
    CGAL_SET_LOG_LEVEL(CGAL::Log::SOME);

    std::vector<S> ns;
    for (unsigned int i=0; i< spheres_.size(); ++i) {
      ns.push_back(S(P(spheres_[i].center().x(),
		       spheres_[i].center().y(),
		       spheres_[i].center().z()),
		     spheres_[i].squared_radius()));
      //std::cout << ns.back() << std::endl;
    }

    typename A::Traits tr(ns.begin(), ns.end());
    A arr(tr);
    typedef typename A::Cross_section CS;
    CS cs;
    
    arr.initialize_at(z_,cs);

    typedef typename CGAL_AOS3_INTERNAL_NS::
      Cross_section_qt_viewer CGAL_AOS3_TARG CSV;
    CSV csv(tr, cs, CGAL::Layer(0));

    typedef typename CGAL_AOS3_INTERNAL_NS::
      Cross_section_qt_viewer_events CGAL_AOS3_TARG CSVE;
    CSVE csve(tr, cs, CGAL::Layer(1));


    typedef typename CGAL_AOS3_INTERNAL_NS::Rational_cross_section CGAL_AOS3_TARG RCS;
    RCS rcs(cs, tr);

    q->set_viewport(up(tr.bbox_3()));
  
    //*q << CGAL::Layer(0);
    //csv(zcur, q);
    //csve(zcur, q);
    
    //q->show_everything();
    //q->redraw();
    int state=0;
    bool finishing=false;
    bool verbose=false;
    bool degen=false;
    while (true) {
      
      csv(zcur, q);
      csve(zcur, q);

      if (state==0) {
	if (!degen) {
	  rcs.set_z(zcur);
	  std::cout << "Auditing arrangment ..." << std::flush;
	  rcs.audit();
	  std::cout << "done" << std::endl;
	} else {
	  degen=false;
	}
      }
      

      std::cout << "Time is " << zcur << std::endl;
      if (state==0) {
	std::cout << "Press return to advance to event: " << std::flush;
      } else if (state==1){
	std::cout << "Press return to process event ";
	if (!cs.visitor().simulator()->empty()) {
	  std::cout << cs.visitor().simulator()->next_event();
	} else {
	  std::cout<<  "None";
	} 
	std::cout << ": "  << std::flush;
      } else {
	std::cout << "Press return to move between events: " << std::flush;
      }
      if (!finishing) {
	char buf[1000];
	std::cin.getline(buf, 1000);
	
	if (buf[0]== 'd') {
	  std::cout << "Time is " << zcur << std::endl;
	  dump(tr, cs);
	} else if (buf[0] == 'f') {
	  finishing=true;
	} else if (buf[0] == 'q') {
	  break;
	} else if (buf[0] == 'v') {
	  verbose=!verbose;
	  if (verbose) {
	    CGAL_SET_LOG_LEVEL(CGAL::Log::LOTS);
	  } else {
	    CGAL_SET_LOG_LEVEL(CGAL::Log::SOME);
	  }
	}
      } else {
	std::cout << std::endl;
      }
      if (state==0) {
	zcur= CGAL::to_interval(cs.visitor().simulator()->next_event_time()).first;
	if (verbose) {
	  std::cout << "Simulator: "  << *cs.visitor().simulator();
	} else {
	  std::cout << "Simulator: "  << cs.visitor().simulator()->current_time();
	  std::cout << std::endl << cs.visitor().simulator()->next_event_time();
	  
	  if (cs.visitor().simulator()->next_event() != Event_key()) {
	    std::cout << ": " << cs.visitor().simulator()->next_event() 
		      << std::endl;
	  }
	}
      } else if (state==1) {
	cs.visitor().simulator()->set_current_event_number(cs.visitor().simulator()->current_event_number()+1 );
	zcur= CGAL::to_interval(cs.visitor().simulator()->current_time()).second;
	if (cs.visitor().simulator()->empty()) {
	  finishing=false;
	} else {
	  if (cs.visitor().simulator()->current_time()
	      == cs.visitor().simulator()->next_event_time()) degen=true;
	}
      } else {
	zcur= .5*(CGAL::to_interval(cs.visitor().simulator()->current_time()).second
		  + CGAL::to_interval(cs.visitor().simulator()->next_event_time()).first);
      }
      state= (state+1)%3;
      
      
    }

    csv(zcur, q);
    csve(zcur, q);
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

  //char buf[1000];
  //char *fname=buf;
  double z;
  if (argc ==1) {
    std::cout << "Enter sweep coordinate: " << std::flush;
    std::cin >> z;
  } else {
    z=std::atof(argv[1]);
  }
  std::string fname;
  if (argc != 3) {
    std::cout << "\nEnter file name: " << std::flush;
    std::cin >> fname;
  } else {
    fname= argv[2];
  }
  std::ifstream in(fname.c_str());
  CGAL_SET_LOG_LEVEL(CGAL::Log::SOME);
  CGAL_AOS3_INTERNAL_NS::read_spheres<Arrangement::Traits::Geom_traits>(in, spheres);
  
 

  typedef Do_work<Arrangement> CS;
  CS cs(spheres, z);
  
  CGAL::Qt_multithreaded_examiner_viewer_2<CS> qtv(cs, argc, argv);
  CGAL::set_error_handler(my_failure_handler);

  return qtv();
}
