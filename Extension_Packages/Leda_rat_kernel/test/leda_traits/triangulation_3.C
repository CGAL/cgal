// provide 3d kernel traits ...
#define CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3

#include <CGAL/basic.h>
#include <CGAL/Triangulation_3.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <list>
#include <iostream> 

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::leda_rat_kernel_traits                           K;
typedef K::Point_3                                             Point;

typedef CGAL::Triangulation_3<K>                               Triang_3;
typedef Triang_3::Vertex_handle                                Vertex_handle;

Triang_3 T;

int main()
{
  std::list<Point>          input;
  std::list<Vertex_handle>  v_handles;
    
  // insert some points
  
  int i;
  Point p;
  
  for(i=0;i<10000;i++){
     p = random_d3_rat_point_in_cube(5000);
     input.push_back(p);
  }
  
  std::list<Point>::const_iterator cit = input.begin();
  for(;cit != input.end(); cit++){
    Vertex_handle vh = T.insert(*cit);
    v_handles.push_back(vh);
  }
  
  
  bool valid = T.is_valid();  
  
  if (valid){
     std::cout << "The 3d triangulation is valid !\n";
  }
  else {
     std::cout << "Error - the 3d triangulation is NOT valid !\n";  
  }
  
  return 0;
}
