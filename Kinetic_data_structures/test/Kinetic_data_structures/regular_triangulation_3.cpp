#define CGAL_CHECK_EXACTNESS
#define CGAL_CHECK_EXPENSIVE
#define CGAL_KINETIC_ENABLE_AUDITING

#include <CGAL/Kinetic/basic.h>

#include <CGAL/Kinetic/Regular_triangulation_3.h>
#include <CGAL/Kinetic/Regular_triangulation_exact_simulation_traits.h>
#include <CGAL/Kinetic/Regular_triangulation_inexact_simulation_traits.h>
#include <CGAL/Kinetic/Regular_triangulation_event_log_visitor_3.h>
#include <algorithm>
#include <iterator>

int main(int argc, char *argv[])
{

 
  if (0) {
    typedef CGAL::Kinetic::Regular_triangulation_inexact_simulation_traits Simulation_traits;
    typedef CGAL::Kinetic::Regular_triangulation_event_log_visitor_3 Visitor;
    typedef CGAL::Kinetic::Regular_triangulation_3<Simulation_traits, Visitor> KDel;
    Simulation_traits simtr(1,1000);
    Simulation_traits::Simulator::Handle sp= simtr.simulator_handle();
      
    KDel kdel(simtr);
    CGAL_SET_LOG_LEVEL(CGAL::Log::NONE);
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
      std::cout << p << std::endl;
      simtr.active_points_3_table_handle()->insert(p); // here 
      ++nread;
      kdel.audit();
    }
    kdel.set_has_certificates(true);
      
    while (sp->next_event_time()
	   < sp->end_time()) {
      sp->set_current_event_number(sp->current_event_number()+1);
    }
  } else {
    typedef CGAL::Kinetic::Regular_triangulation_exact_simulation_traits Simulation_traits;
    typedef CGAL::Kinetic::Regular_triangulation_event_log_visitor_3 Visitor;
    typedef CGAL::Kinetic::Regular_triangulation_3<Simulation_traits, Visitor> KDel;
    typedef KDel::Triangulation IT;
    Simulation_traits simtr(1,1000);
    Simulation_traits::Simulator::Handle sp= simtr.simulator_handle();

    typedef CGAL::Regular_triangulation_3<Simulation_traits::Static_kernel> ST;
    typedef ST::Weighted_point SWP;
    typedef ST::Bare_point BP;

    Simulation_traits isimtr(1,1000);
    IT it(isimtr.instantaneous_kernel_object());
    it.geom_traits().set_time(Simulation_traits::Instantaneous_kernel::NT(0));
    
    ST st;
    KDel kdel(simtr);
    kdel.triangulation().geom_traits().set_time(Simulation_traits::Instantaneous_kernel::NT(0));
    if (argc==1) {
      CGAL_SET_LOG_LEVEL(CGAL::Log::NONE);
      std::ifstream in("data/regular_triangulation_3.input");
      if (!in) {
	std::cerr << "Error opening input file: " << "data/regular_triangulation_3.input" 
		  << std::endl;
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
    } else {
      CGAL_SET_LOG_LEVEL(CGAL::Log::LOTS);
       std::ifstream in(argv[1]);
       if (!in) {
	 std::cerr << "Error opening input file: " <<argv[1] << std::endl;
	 return EXIT_FAILURE;
       }
       char buf[1000];
       int nread=0;
       std::map<SWP, int> pm;
       pm[st.vertices_begin()->point()]= 0;
       while (true ) {
	 in.getline(buf, 1000);
	 if (!in) break;
	 std::istringstream il(buf);
	 Simulation_traits::Kinetic_kernel::Weighted_point_3 p;
	 il >> p;
	 if (0) {
	   Simulation_traits::Instantaneous_kernel::Point_3 ip
	     =isimtr.active_points_3_table_handle()->insert(p);
	   IT::Cell_handle h= it.locate(ip);
	   if (h!= IT::Cell_handle() 
	       && h->vertex(0) != IT::Vertex_handle()
	       && h->vertex(1) != IT::Vertex_handle()
	       && h->vertex(2) != IT::Vertex_handle()
	       && h->vertex(3) != IT::Vertex_handle()) {
	     std::cout << "Instant located in " << h->vertex(0)->point() << " "
		       << h->vertex(1)->point() << " "
		       << h->vertex(2)->point() << " "
		       << h->vertex(3)->point() << std::endl;
	   } else {
	     std::cout << "Instant located outside hull" <<std::endl; 
	   }
	   std::vector<IT::Facet> bfacets;
	   std::vector<IT::Cell_handle> cells;
	   std::vector<IT::Facet> ifacets;
	   if (it.dimension() == 3) {
	     it.find_conflicts(ip, h, back_inserter(bfacets), 
			       back_inserter(cells),back_inserter(ifacets));
	   }
	   
	   it.insert(ip, h);
	 }
	 if (0) {
	   Simulation_traits::Static_kernel::FT ct(0);
	   SWP sp(BP(p.point().x()(ct), p.point().y()(ct),
		     p.point().z()(ct)), p.weight()(ct));
	   ST::Cell_handle h= st.locate(sp);
	   if (h!= ST::Cell_handle()
	       && h->vertex(0) != ST::Vertex_handle()
	       && h->vertex(1) != ST::Vertex_handle()
	       && h->vertex(2) != ST::Vertex_handle()
	       && h->vertex(3) != ST::Vertex_handle()) {
	     std::cout << "Static located in " << pm[h->vertex(0)->point()]-1 << " "
		       << pm[h->vertex(1)->point()]-1 << " "
		       << pm[h->vertex(2)->point()]-1 << " "
		       << pm[h->vertex(3)->point()]-1 << std::endl;
	   } else {
	     std::cout << "Instant located outside hull" <<std::endl; 
	   }
	   st.insert(sp, h);
	   pm[sp]= pm.size()-1;
	 }
	 simtr.active_points_3_table_handle()->set_is_editing(true);
	 Simulation_traits::Active_points_3_table::Key vp=simtr.active_points_3_table_handle()->insert(p);
	 {
	   
	   KDel::Cell_handle h=  kdel.triangulation().locate(vp);
	   if (h!= KDel::Cell_handle() 
	       && h->vertex(0) != KDel::Vertex_handle()
	       && h->vertex(1) != KDel::Vertex_handle()
	       && h->vertex(2) != KDel::Vertex_handle()
	       && h->vertex(3) != KDel::Vertex_handle()) {
	     std::cout << "Kinetic XT located in " << h->vertex(0)->point() << " "
		       << h->vertex(1)->point() << " "
		       << h->vertex(2)->point() << " "
		       << h->vertex(3)->point() << std::endl;
	   } else {
	     std::cout << "Instant located outside hull" <<std::endl; 
	   }
	 }
	 simtr.active_points_3_table_handle()->set_is_editing(false);
	 //std::cout << p << std::endl;
	  // here 
	 ++nread;
	 kdel.audit();
       }
    }
    kdel.set_has_certificates(true);
    if (simtr.simulator_handle()->has_audit_time()) kdel.audit();

    while (sp->next_event_time()
	   < sp->end_time()) {
      sp->set_current_event_number(sp->current_event_number()+1);
    }

    if (argc==1) {
      std::ifstream out("data/regular_triangulation_3.output");
      if (!out) {
	std::cerr << "Error opening input file: " << "data/regular_triangulation_3.output" << std::endl;
	return EXIT_FAILURE;
      }
      
      std::vector<std::string> expected;
      while (true) {
	char buf[1000];
	out.getline(buf, 1000);
	if (!out) break;
	expected.push_back(buf);
      }
      
      int error_count=0;
      unsigned int line=0;
      for (CGAL::Kinetic::Delaunay_triangulation_event_log_visitor_3::Event_iterator it = kdel.visitor().events_begin();
	   it != kdel.visitor().events_end(); ++it) {
	if (expected.size() <= line || *it != expected[line]) {
	  //std::cerr << "ERROR Got event: " << *it << std::endl;
	  //std::cerr << "      Expected event: " << buf << std::endl;
	  ++error_count;
	}
	++line;
      }
      if (error_count != 0) {


	std::copy(kdel.visitor().events_begin(), kdel.visitor().events_end(),
		  std::ostream_iterator<std::string>(std::cout, "\n"));

	line=0;
	error_count=0;
	for (CGAL::Kinetic::Delaunay_triangulation_event_log_visitor_3::Event_iterator it = kdel.visitor().events_begin();
	     it != kdel.visitor().events_end(); ++it) {
	  if (expected.size() <= line || *it != expected[line]) {
	    int found_line=-1;
	    for (unsigned int i=line; i < expected.size(); ++i) {
	      if (expected[i] == *it) {
		found_line=i;
		break;
	      }
	    }
	    if (found_line != -1) {
	      for (int i=line; i< found_line; ++i) {
		std::cerr << "ERROR Missing event: " << expected[i] << std::endl;
		++error_count;
	      }
	      assert(*it == expected[found_line]);
	      std::cerr << "Matched Event: " << *it << std::endl;
	      line=found_line+1; 
	    } else {
	      std::cerr << "ERROR Extra event: " << *it << std::endl;
	      ++error_count;
	    }
	  } else {
	    std::cerr << "Matched event: " << *it << std::endl;
	    ++line;
	  }
	}
	for (; line != expected.size(); ++line) {
	  std::cerr << "ERROR Missing event: " << expected[line] << std::endl;
	  ++error_count;
	}
      }
      if (error_count != 0) {
	std::cerr << "ERROR " << error_count << " errors in " << kdel.visitor().size() << " events.\n";
      }
    }
  

   
    CGAL::Kinetic::internal::write_debug_counters(std::cout);
    if (CGAL::Kinetic::internal::get_static_audit_failures() != 0) return EXIT_FAILURE;
    else return EXIT_SUCCESS;
  }
}
