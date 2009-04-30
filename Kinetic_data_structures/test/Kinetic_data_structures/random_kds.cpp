#define CGAL_CHECK_EXPENSIVE
#define CGAL_CHECK_EXACTNESS
#define CGAL_KINETIC_DO_NOT_USE_LAZY_EXACT

#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/Regular_triangulation_exact_simulation_traits.h>
#include <CGAL/Kinetic/Sort.h>
#include <CGAL/Kinetic/Delaunay_triangulation_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_3.h>
#include <CGAL/Kinetic/Regular_triangulation_3.h>
#include <CGAL/Random.h>

template <class P>
P rp(CGAL::Random &rand, int deg) {
  std::vector<typename P::NT> coefs(deg+1);
  for (int i=0; i< deg+1; ++i){
    double mag= 10.0/(i+1);
    coefs[i]= rand.get_double(-mag, mag);
  }
  return P(coefs.begin(), coefs.end());
}

int main(int, char *[]) {

  CGAL::Random r(time(NULL));

  int which = r.get_int(0,4);
  double nump= r.get_int(25,50);
  double end_time= r.get_double(10,100);
  if (which ==0) {
    std::cout << "Checking Delaunay_2" << std::endl;
    typedef CGAL::Kinetic::Exact_simulation_traits Traits;
    typedef CGAL::Kinetic::Delaunay_triangulation_2<Traits> DS_;
    typedef Traits::Kinetic_kernel::Motion_function F;
    typedef Traits::Kinetic_kernel::Point_2 Point;
    Traits tr(0, end_time);
    DS_ ds(tr);
    std::cout << "Points are:\n";

    for (int i=0; i< nump/1.5; ++i){
      Point pt(rp<F>(r, r.get_int(0,5)),
	       rp<F>(r, r.get_int(0,5)));
      std::cout << pt << std::endl;
      tr.active_points_2_table_handle()->insert(pt);
       if (i == 10) {
	ds.set_has_certificates(true);
	ds.audit();
	ds.set_has_certificates(false);
	tr.simulator_handle()->audit_events();
      }
    }
    ds.set_has_certificates(true);
    ds.audit();
    tr.simulator_handle()->set_current_time(tr.simulator_handle()->end_time());
    ds.audit();
    std::cout << "Processed " << tr.simulator_handle()->current_event_number() << " events"<<std::endl;
  } else if (which ==1) {
    std::cout << "Checking Sort" << std::endl;

    typedef CGAL::Kinetic::Exact_simulation_traits Traits;
    typedef CGAL::Kinetic::Sort<Traits> DS_;
    typedef Traits::Kinetic_kernel::Motion_function F;
    typedef Traits::Kinetic_kernel::Point_1 Point;
    Traits tr(0, end_time);
    DS_ ds(tr);
    std::cout << "Points are:\n";
    for (int i=0; i< nump; ++i){
      Point pt(rp<F>(r, r.get_int(0,5)));
      std::cout << pt << std::endl;
      tr.active_points_1_table_handle()->insert(pt);
      if (i == 10) {
	//ds.set_has_certificates(true);
	ds.audit();
	//ds.set_has_certificates(false);
	tr.simulator_handle()->audit_events();
      }
    }
    std::cout << *tr.active_points_1_table_handle() << std::endl;
    ds.audit();
    //ds.set_has_certificates(true);
    tr.simulator_handle()->set_current_time(tr.simulator_handle()->end_time());
    ds.audit();
    std::cout << "Processed " << tr.simulator_handle()->current_event_number() << " events"<<std::endl;

  } else if (which == 2) {
    std::cout << "Checking Delaunay_3" << std::endl;
    typedef CGAL::Kinetic::Exact_simulation_traits Traits;
    typedef CGAL::Kinetic::Delaunay_triangulation_3<Traits> DS_;
    typedef Traits::Kinetic_kernel::Motion_function F;
    typedef Traits::Kinetic_kernel::Point_3 Point;
    Traits tr(0, end_time);
    DS_ ds(tr);
    std::cout << "Points are:\n";
    for (int i=0; i< nump/2; ++i){
      Point pt(rp<F>(r, r.get_int(0,3)),
	       rp<F>(r, r.get_int(0,3)),
	       rp<F>(r, r.get_int(0,3)));
      std::cout << pt << std::endl;
      tr.active_points_3_table_handle()->insert(pt);
      if (i == 10) {
	ds.set_has_certificates(true);
	ds.audit();
	ds.set_has_certificates(false);
	tr.simulator_handle()->audit_events();
      }
    }
    std::cout << *tr.active_points_3_table_handle() << std::endl;
    ds.set_has_certificates(true);
    ds.audit();
    tr.simulator_handle()->set_current_time(tr.simulator_handle()->end_time());
    ds.audit();
    std::cout << "Processed " << tr.simulator_handle()->current_event_number() << " events"<<std::endl;
    
  } else {
    std::cout << "Checking regular_3" << std::endl;
    typedef CGAL::Kinetic::Regular_triangulation_exact_simulation_traits Traits;
    typedef CGAL::Kinetic::Regular_triangulation_3<Traits> DS_;
    typedef Traits::Kinetic_kernel::Motion_function F;
    typedef Traits::Kinetic_kernel::Point_3 Bare_point;
    typedef Traits::Kinetic_kernel::Weighted_point_3 Point;
    Traits tr(0, end_time);
    DS_ ds(tr);
    std::cout << "Points are:\n";

    for (int i=0; i< nump/2; ++i){
      Point pt = Point(Bare_point(rp<F>(r, r.get_int(0,5)),
			          rp<F>(r, r.get_int(0,5)),
			          rp<F>(r, r.get_int(0,5))),
	               rp<F>(r, r.get_int(0,5)));
      std::cout << pt << std::endl;
      tr.active_points_3_table_handle()->insert(pt);
      ds.audit();
      if (i == 10) {
	ds.set_has_certificates(true);
	ds.audit();
	ds.set_has_certificates(false);
	tr.simulator_handle()->audit_events();
      }
    }
    std::cout << *tr.active_points_3_table_handle() << std::endl;
    ds.set_has_certificates(true);
    ds.audit();
    tr.simulator_handle()->set_current_time(tr.simulator_handle()->end_time());
    ds.audit();
    std::cout << "Processed " << tr.simulator_handle()->current_event_number() << " events"<<std::endl;
  }
  
  return EXIT_SUCCESS;
}
