//#define CGAL_NO_STATIC_FILTERS

#include <CGAL/basic.h>

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Gmpq.h>

//robust classes
#include <CGAL/Regular_triangulation_sphere_traits_2.h>
#include <CGAL/Regular_triangulation_on_sphere_2.h>
//#include <CGAL/Triangulation_sphere_traits_2.h>
#include <CGAL/Triangulation_on_sphere_2.h>

//manuel prototypes
//#include <CGAL/Regular_triangulation_sphere_traits_2_manuel_proto.h>
//#include <CGAL/Regular_triangulation_on_sphere_2_manuel_proto.h>
//#include <CGAL/Triangulation_on_sphere_2_manuel_proto.h>


#include <CGAL/Random.h>
#include <CGAL/Timer.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Regular_triangulation_sphere_traits_2< K > FK;
template class CGAL::Regular_triangulation_on_sphere_2<FK>;
typedef CGAL::Regular_triangulation_on_sphere_2<FK>  RTOS;

//typedef CGAL::Regular_triangulation_sphere_traits_2_manuel_proto< K > FK2;
//template class CGAL::Regular_triangulation_on_sphere_2_manuel_proto<FK2>;
//typedef CGAL::Regular_triangulation_on_sphere_2_manuel_proto<FK2> RTOSPROTO;


void compute_times()
{	
	int no_of_pts;

	typedef K::Point_3 Point_3;
  std::list<Point_3> pts3;

 // std::cin >> no_of_pts;
	no_of_pts = 1000;
  for (int count=0; count<no_of_pts; count++) {
    Point_3 p;
		std::cin >> p;
    if(p != Point_3(0,0,0)) pts3.push_back(p);
		else count--;
  }
  
  std::cout << "REGULAR TRIANGULATION ON SPHERE" << std::endl; 
  {
 	  RTOS rtos;
	  //rtos.insert_four_init_vertices();
    {
      CGAL::Timer timer_reg; timer_reg.start();
      rtos.insert(pts3.begin(),pts3.end());
      timer_reg.stop();
			std::cout << timer_reg.time() << std::endl;
    }
		if(!rtos.is_valid()) { std::cout << "problem!" << std::endl;  }
  }
  
  
  std::cout << "REGULAR TRIANGULATION ON SPHERE PROTO" << std::endl; 
  {
	  RTOS rtosprot;
    {
      CGAL::Timer timer_reg; timer_reg.start();
      rtosprot.insert(pts3.begin(),pts3.end());
      timer_reg.stop();
			std::cout << timer_reg.time() << std::endl;
    }
		if(!rtosprot.is_valid()) { std::cout << "problem!" << std::endl; exit(0); }
  }

	 
}

int main(int argc, char* argv[]) 
{
	compute_times();
	return 0;
}

