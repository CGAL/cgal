#include <CGAL/KDS/basic.h>

#include <CGAL/KDS/Delaunay_triangulation_2.h>
#include <CGAL/KDS/Insert_event.h>
#include <CGAL/Random.h>
#include <CGAL/KDS/Exact_simulation_traits_2.h>

int main(int, char *[]){

  typedef CGAL::KDS::Exact_simulation_traits_2 Simulation_traits;
  typedef Simulation_traits::Kinetic_kernel::Point_2 Moving_point_2;

  typedef CGAL::KDS::Delaunay_triangulation_2<Simulation_traits> KDel;
  
  Simulation_traits simtr;
  Simulation_traits::Simulator::Pointer sp= simtr.simulator_pointer();

  KDel kdel(simtr);

  CGAL::Random rand(std::time(NULL));
  
  Simulation_traits::Function_kernel::Construct_function cf
    = simtr.function_kernel_object().construct_function_object();
  
  for (unsigned int i=0; i< 8; ++i){
    Moving_point_2 mp(cf(rand.get_double()*10-5, 
			 rand.get_double(),
			 rand.get_double()),
	  cf(rand.get_double()*10-5, 
	     rand.get_double(), 
	     rand.get_double()));
    Simulation_traits::Simulator::Time t(static_cast<double>(i)/4.0);
    
    simtr.moving_point_table_pointer()->insert(mp);
  }

  while (simtr.simulator_pointer()->next_event_time() 
	 < simtr.simulator_pointer()->end_time()){
    sp->set_current_event_number(sp->current_event_number()+1);
    std::cout << "At time " << simtr.simulator_pointer()->current_time() << ":\n";
    std::cout << kdel.triangulation_data_structure();
  }

  return EXIT_SUCCESS;
};
