//example6

#define NAIVE_POINT_LOCATION
#define CGAL_NO_PM_DEFAULT_POINT_LOCATION

//#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Pm_segment_exact_traits.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

#include <CGAL/IO/Planar_map_iostream.h>

#ifdef NAIVE_POINT_LOCATION
#include <CGAL/Pm_naive_point_location.h>
#endif

//typedef CGAL::Homogeneous<long>                    coord_t;
typedef CGAL::Cartesian<double>                      Coord_t;
typedef CGAL::Pm_segment_exact_traits<Coord_t>       pmtraits;
typedef pmtraits::Point                              point;
typedef pmtraits::X_curve                            curve;
typedef CGAL::Pm_default_dcel<pmtraits>              pmdcel;
typedef CGAL::Planar_map_2<pmdcel,pmtraits>          Pmap;

int main()
{
  // creating an instance of CGAL::Planar_map_2<pmdcel,pmtraits>
#ifndef NAIVE_POINT_LOCATION   
  Pmap pm;
#else
  CGAL::Pm_naive_point_location< Pmap > naive_pl;
  Pmap pm(&naive_pl);
#endif
	
  curve cv[6];
  int i;
	
  CGAL::set_ascii_mode(std::cout);
	
  point a1(1, 1), a2(1, 0), a3(0, 0), a4(0, 1), a5(1,4,2) ;
	
  /*
	a5 
	/\  
	/  \
	a4  ----  a1
	|    | 
	|    | 
	|    |
	a3  ----  a2
	*/        
	
	
	// those curves are about to enter to pm
  cv[0] = curve(a1, a2);
  cv[1] = curve(a2, a3);
  cv[2] = curve(a3, a4);
  cv[3] = curve(a4, a5);
  cv[4] = curve(a5, a1);
  cv[5] = curve(a1, a4);
	
  std::cout << "the curves of the map :" << std::endl; 
  for (i = 0; i < 6; i++)
    std::cout << cv[i] << std::endl;
	
  std::cout << std::endl;
	
	
  // insert the five curves to the map
  std::cout << "inserting the curves to the map..." << std::endl;
  CGAL::Planar_map_2<pmdcel,pmtraits>::Halfedge_handle e[6];  
	
  e[0]=pm.insert(cv[0]);
	
  for (i = 1; i < 4; i++)
    {
      e[i]=pm.insert_from_vertex(cv[i],e[i-1]->target(), true);
    }
	
  e[4]=pm.insert_at_vertices(cv[4],e[0]->source(),e[3]->target() );
	
  e[5]=pm.insert_at_vertices(cv[5],e[0]->source(),e[2]->target() );
	
  /*             
	e3  /\  e4
	/  \
	----
	| e5 | 
	e2 |    | e0
	|    |
	----
	e1 
	*/
	
	// check the validity of the map
  std::cout << "check map validity... ";
  if (pm.is_valid())
    std::cout << "map valid!" << std::endl;
  else
    std::cout << "map invalid!" << std::endl;
  std::cout << std::endl << "Printing map: " << std::endl;
	
  std::cout << pm ;
	
  return 0;  
}
