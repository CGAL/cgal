//example2

#include <CGAL/Homogeneous.h>             //change from example1
#include <CGAL/Pm_segment_exact_traits.h> //change from example1

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

using namespace CGAL;

typedef Homogeneous<long>                 coord_t;  //change from example1
typedef Pm_segment_exact_traits<coord_t>  pmtraits; //change from example1

typedef pmtraits::Point                   point;
typedef pmtraits::X_curve                 curve;
typedef Pm_default_dcel<pmtraits>         pmdcel;

int main()
{
  // creating an instance of Planar_map_2<pmdcel,pmtraits>
  //Pm_naive_point_location_strategy<pmdcel,pmtraits> pl_strategy;  
  //Planar_map_2<pmdcel,pmtraits> pm(&pl_strategy);
  Planar_map_2<pmdcel,pmtraits> pm;

  curve cv[5];
  int i;

  point a1(100, 0), a2(20, 50), a3(180, 50), a4(100, 100);

  // those curves are about to be inserted to pm
  cv[0] = curve(a1, a2);
  cv[1] = curve(a1, a3);
  cv[2] = curve(a2, a3);
  cv[3] = curve(a2, a4);
  cv[4] = curve(a3, a4);
  
  std::cout << "the curves of the map :" << std::endl; 
  for (i = 0; i < 5; i++)
    std::cout << cv[i] << std::endl;

  std::cout << std::endl;

  // insert the five curves to the map
  std::cout << "inserting the curves to the map..." << std::endl;
  for (i = 0; i < 5; i++)
  {
    std::cout << "inserting curve" << i << std::endl;
    pm.insert(cv[i]);
  }
  
  // check the validity of the map
  std::cout << "check map validity... ";
  if (pm.is_valid())
    std::cout << "map valid!" << std::endl;
  else
    std::cout << "map invalid!" << std::endl;
  std::cout << std::endl;
  
  // vertical ray shooting upward from p
  point p(95, 30);
  Planar_map_2<pmdcel,pmtraits>::Halfedge_handle e;  
  Planar_map_2<pmdcel,pmtraits>::Locate_type lt;

  std::cout << std::endl << "upward vertical ray shooting from " << p;
  std::cout << std::endl; 

  e=pm.vertical_ray_shoot(p, lt, true);
  std::cout << "returned the curve :" << e->curve() << " oriented toward " 
            << e->target()->point() << std::endl;

  //testing the removal function

  std::cout << "\nremoving the edge ... " << std::endl;

  pm.remove_edge(e);

  if (pm.is_valid())
    std::cout << "map valid!" << std::endl;
  else
    std::cout << "map invalid!" << std::endl;
  std::cout << std::endl;


  return 0;  
}

