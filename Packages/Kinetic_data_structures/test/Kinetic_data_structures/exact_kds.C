#include <CGAL/KDS/Sort.h>
#include <CGAL/KDS/Insert_event.h>
#include <CGAL/KDS/Exact_simulation_traits_1.h>
#include <cstdlib>
#include "include/sort_test.h"



int main(int argc, char *argv[]){
  unsigned int num_points=20;
  unsigned int degree =2;
  if (argc==3){
    num_points= std::atoi(argv[1]);
    degree= std::atoi(argv[2]);
  }
  std::cout << "Using " << num_points  << " degree " << degree << " points.\n";
  //CGAL_KDS_SET_LOG_LEVEL(CGAL::KDS::LOG_SOME);
  typedef CGAL::KDS::Exact_simulation_traits_1 Tr;
  Tr tr;
  typedef Tr::Simulator::Time Time;

  typedef CGAL::KDS::Insert_event<Tr::Moving_point_table> MOI;
  typedef Tr::Kinetic_kernel::Point_2 MP;
  
  for (unsigned int i=0; i< num_points; ++i){
    std::vector<double> coefs;
    for (unsigned int j=0; j<= degree; ++j){
      coefs.push_back(static_cast<double>(std::rand())/static_cast<double>(RAND_MAX));
    }
    tr.simulator_pointer()->new_event(Time(i/100.0), MOI(MP(Tr::Function_kernel::Function(coefs.begin(), coefs.end()), 
							    Tr::Function_kernel::Function(0)), tr.moving_point_table_pointer()));
  }

  

  bool error= sort_test<Tr>(tr);
  if (error){
    return EXIT_FAILURE;
  } else {
    return EXIT_SUCCESS;
  }
};
