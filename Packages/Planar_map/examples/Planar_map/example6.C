// examples/Planar_map/example6.C
// ------------------------------
#define NAIVE_POINT_LOCATION
#define CGAL_NO_PM_DEFAULT_POINT_LOCATION

//#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Pm_segment_exact_traits.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

//#include <CGAL/IO/Planar_map_iostream.h>
#include <CGAL/IO/Pm_file_writer.h>
#include <CGAL/IO/write_pm.h>

#ifdef NAIVE_POINT_LOCATION
#include <CGAL/Pm_naive_point_location.h>
#endif

using CGAL::write_pm; 

//typedef CGAL::Homogeneous<long>                    coord_t;
typedef CGAL::Cartesian<double>                      Coord_t;
typedef CGAL::Pm_segment_exact_traits<Coord_t>       Pmtraits;
typedef Pmtraits::Point                              Point;
typedef Pmtraits::X_curve                            Curve;
typedef CGAL::Pm_default_dcel<Pmtraits>              Pmdcel;
typedef CGAL::Planar_map_2<Pmdcel, Pmtraits>         Planar_map;
typedef CGAL::Pm_file_writer<Planar_map>             Pm_writer;

int main()
{
  // creating an instance of CGAL::Planar_map_2<pmdcel,pmtraits>
#ifndef NAIVE_POINT_LOCATION   
  Planar_map pm;
#else
  CGAL::Pm_naive_point_location< Planar_map > naive_pl;
  Planar_map pm(&naive_pl);
#endif
  Pm_writer verbose_writer(std::cout, pm, true);
	
  Curve cv[6];
  int i;
	
  CGAL::set_ascii_mode(std::cout);
	
  Point a1(1, 1), a2(1, 0), a3(0, 0), a4(0, 1), a5(1,4,2) ;
	
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
  cv[0] = Curve(a1, a2);
  cv[1] = Curve(a2, a3);
  cv[2] = Curve(a3, a4);
  cv[3] = Curve(a4, a5);
  cv[4] = Curve(a5, a1);
  cv[5] = Curve(a1, a4);
	
  std::cout << "the curves of the map :" << std::endl; 
  for (i = 0; i < 6; i++)
    std::cout << cv[i] << std::endl;
	
  std::cout << std::endl;
	
	
  // insert the five curves to the map
  std::cout << "inserting the curves to the map..." << std::endl;
  Planar_map::Halfedge_handle e[6];  
	
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
  std::cout << std::endl << "* * * Printing map: " << std::endl << std::endl;
	
  write_pm(pm, verbose_writer, std::cout);
       
  return 0;  
}
