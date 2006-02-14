
#include <CGAL/basic.h>
#include <CGAL/Triangulation_2.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <LEDA/random_rat_point.h>

#include <list>
#include <iostream> 

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::leda_rat_kernel_traits                           K;
typedef K::Point_2                                             Point;

typedef CGAL::Triangulation_2<K>                               Triang_2;
typedef Triang_2::Vertex_handle                                Vertex_handle;

Triang_2 T;

int main()
{
  std::list<Point>          input;
  std::list<Vertex_handle>  v_handles;
  
  // insert some points
  
  int i;
  Point p;
  
  for(i=0;i<10000;i++){
     random_point_in_square(p, 5000);
     input.push_back(p);
  }
  
  std::list<Point>::const_iterator cit = input.begin();
  for(;cit != input.end(); cit++){
    Vertex_handle vh = T.insert(*cit);
    v_handles.push_back(vh);
  }
  
  // delete some points
 
  std::list<Vertex_handle>::const_iterator vit = v_handles.begin();  
  for(i=0;i<100;i++){
     T.remove(*vit);
     vit++;
  }  
  
  bool valid = T.is_valid();
  
  if (valid){
     std::cout << "The triangulation is valid !\n";
  }
  else {
     std::cout << "Error - the triangulation is NOT valid !\n";  
  }
  
  return 0;
}
