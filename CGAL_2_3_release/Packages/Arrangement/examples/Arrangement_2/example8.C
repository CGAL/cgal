//examples/Arrangement_2/example8.C

// Define shorter names to please linker (g++/egcs)
#define Arrangement_2 Ar
#define _In_place_list_iterator IPLI
#define Homogeneous Ho
#define Arr_segment_exact_traits ASET

#include <CGAL/Homogeneous.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>
//#include <CGAL/leda_integer.h>

#include <CGAL/Arr_segment_exact_traits.h>

//typedef leda_integer                                NT;
typedef long                                        NT;

typedef CGAL::Homogeneous<NT>                       R;

typedef CGAL::Arr_segment_exact_traits<R>           Traits;

typedef Traits::Point                               Point;
typedef Traits::X_curve                             X_curve;

typedef CGAL::Arr_base_node<X_curve>                Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>            Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node>  Arr_2;

int main(int argc, char* argv[])
{
  Arr_2 arr; 

  int num_curves;
  int x,y;
  std::cout << "Enter number of segments: " ;
  std::cin >> num_curves;
  while (num_curves--) {
    std::cout << "Enter source coordinates (2 integers): " ;
    std::cin >> x >> y;
    Point s(x,y);

    std::cout << "Enter target coordinates (2 integers): " ;
    std::cin >> x >> y;
    Point t(x,y);

    X_curve seg(s,t);
    arr.insert(seg);
  }

  std::cout << "Enter point for ray shooting (2 integers): " ;
  std::cin >> x >> y;
  Point p(x,y);

  Arr_2::Halfedge_handle e=arr.halfedges_begin();
  Arr_2::Locate_type lt;
  e = arr.vertical_ray_shoot(p,lt,true);
  
  if (lt==Arr_2::UNBOUNDED_FACE) {
    std::cout << "UNBOUNDED_FACE" << std::endl;
  }
  else {
    std::cout << "The halfedge shot is :\n";
    std::cout << "(Using homogeneous coordinates <hx, hy, hw>, ";
    std::cout << "where <x, y>=<hx/hw, hy/hw>)" << std::endl;
    std::cout << e->source()->point() << " -> " << e->target()->point();
    std::cout << std::endl;
  }

  return 0;  
}
