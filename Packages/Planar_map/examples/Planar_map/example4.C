//example4

#include <CGAL/Homogeneous.h>
#include <CGAL/Pm_segment_exact_traits.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

//for operator<<
#include <CGAL/IO/Planar_map_iostream.h>

using namespace CGAL;

typedef Homogeneous<long>                  coord_t;
typedef Pm_segment_exact_traits<coord_t>   Pmtraits;

typedef Pmtraits::Point                    Point;
typedef Pmtraits::X_curve                  Curve;

typedef Pm_default_dcel<Pmtraits>          Pmdcel;

int main()
{
  // creating an instance of Planar_map_2<Pmdcel,Pmtraits>
  Planar_map_2<Pmdcel,Pmtraits> pm;

  Curve cv[5];
  int i;

  set_ascii_mode(std::cout);

  Point a1(100, 0), a2(20, 50), a3(180, 50), a4(100, 100);

  // those curves are about to be inserted to pm
  cv[0] = Curve(a1, a2);
  cv[1] = Curve(a1, a3);
  cv[2] = Curve(a2, a3);
  cv[3] = Curve(a2, a4);
  cv[4] = Curve(a3, a4);
  
  Planar_map_2<Pmdcel,Pmtraits>::Halfedge_handle e[5];  
  // insert the five curves to the map and return e[i]
  for (i = 0; i < 5; i++)
  {
    e[i]=pm.insert(cv[i]);
    std::cout << "is " ;
    if (!pm.is_valid() )
      std::cout << "in" ;
    std::cout << "valid" << std::endl  ;
  }
  
   //map before splitting the edge and adding curve
  std::cout << "Map before:" << std::endl;
  std::cout << pm << std::endl;
  

  //splitting e[2] of the map at the middle and inserting an edge between the 
  // new vertex and the vertex at a1
  Point p(100, 50);
  Curve c1(a2,p);
  Curve c2(p,a3);

  Planar_map_2<Pmdcel,Pmtraits>::Halfedge_handle se = pm.split_edge(e[2],c1,c2); 

  pm.insert_at_vertices( Curve(p,a1), se->target(),e[0]->source() );

  std::cout << "Map after:" << std::endl;
  std::cout << pm;

  return 0;  
}
