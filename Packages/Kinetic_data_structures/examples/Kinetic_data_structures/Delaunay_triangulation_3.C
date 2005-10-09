#include <CGAL/KDS/Exact_simulation_traits_3.h>
#include <CGAL/KDS/Delaunay_triangulation_3.h>
#include <CGAL/KDS/Delaunay_triangulation_event_log_visitor_3.h>

int main(int , char *[]){

  typedef CGAL::KDS::Exact_simulation_traits_3 Traits;
  typedef CGAL::KDS::Delaunay_triangulation_3<Traits, 
    CGAL::KDS::Delaunay_triangulation_event_log_visitor_3> KDel;

  Traits tr;
  KDel kdel(tr);
  
 
  Traits::Simulator::Pointer sp= tr.simulator_pointer();

  Traits::Function_kernel::Construct_function cf
    = tr.function_kernel_object().construct_function_object();
  
  CGAL::Random rand;
  for (int i=0; i< 10; ++i){
    std::vector<double> coefsx, coefsy, coefsz;
    for (int j=0; j< 4; ++j){
      coefsx.push_back((rand.get_double()*10-5)/(j+1));
      coefsy.push_back((rand.get_double()*10-5)/(j+1));
      coefsz.push_back((rand.get_double()*10-5)/(j+1));
    }
    Traits::Kinetic_kernel::Point_3 mp(Traits::Function_kernel::Function(coefsx.begin(), coefsx.end()),
				       Traits::Function_kernel::Function(coefsy.begin(), coefsy.end()),
				       Traits::Function_kernel::Function(coefsz.begin(), coefsz.end()));
    //std::cout << mp << std::endl;
    tr.moving_point_table_pointer()->insert(mp);
  }

  kdel.set_has_certificates(true);

  sp->set_current_time(sp->end_time());

  std::copy(kdel.visitor().begin(), kdel.visitor().end(),
	    std::ostream_iterator<std::string>(std::cout, "\n"));
  return EXIT_SUCCESS;
};
