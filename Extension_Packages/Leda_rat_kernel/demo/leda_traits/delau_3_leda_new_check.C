// we use the new style kernel checker ...

// provide 3d kernel traits ...
#define CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3
#define CGAL_NO_DEPRECATED_CODE

#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Kernel_checker.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <LEDA/integer.h>
#include <iostream>
#include <list>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif


// LEDA -> CGAL conversion

CGAL_BEGIN_NAMESPACE

struct leda_cgal_converter_traits {
  typedef leda_to_cgal_2   Converter_object_2;
  typedef leda_to_cgal_3   Converter_object_3;          

  typedef leda_to_cgal_2   Converter_point_2;
  typedef leda_to_cgal_2   Converter_vector_2;
  typedef leda_to_cgal_2   Converter_direction_2;
  typedef leda_to_cgal_2   Converter_line_2;
  typedef leda_to_cgal_2   Converter_ray_2;
  typedef leda_to_cgal_2   Converter_segment_2;
  typedef leda_to_cgal_2   Converter_triangle_2;
  typedef leda_to_cgal_2   Converter_circle_2;
  typedef leda_to_cgal_2   Converter_iso_rectangle_2;

  typedef leda_to_cgal_3   Converter_point_3;
  typedef leda_to_cgal_3   Converter_vector_3;
  typedef leda_to_cgal_3   Converter_direction_3;
  typedef leda_to_cgal_3   Converter_line_3;
  typedef leda_to_cgal_3   Converter_plane_3;
  typedef leda_to_cgal_3   Converter_ray_3;
  typedef leda_to_cgal_3   Converter_segment_3;
  typedef leda_to_cgal_3   Converter_triangle_3;
  typedef leda_to_cgal_3   Converter_tetrahedron_3;
  typedef leda_to_cgal_3   Converter_sphere_3;
  typedef leda_to_cgal_3   Converter_iso_cuboid_3; 
};

CGAL_END_NAMESPACE


typedef CGAL::leda_rat_kernel_traits                    K1;
typedef CGAL::Homogeneous<leda_integer>                 K2;

typedef CGAL::leda_cgal_converter_traits                CONV;

typedef CGAL::Kernel_checker<K1, K2, CONV>              K;

typedef CGAL::Delaunay_triangulation_3<K>               DT;
typedef K::Point_3                                      Point;


int main()
{
  DT T;

  // insertion of points
  int x,y,z;
  std::list<Point> L2;
  
  std::cout << "before construction !\n"; std::cout.flush();  
  
  for (z=0 ; z<10 ; z++)
    for (y=0 ; y<10 ; y++)
      for (x=0 ; x<10 ; x++) {
        leda_d3_rat_point pact(x,y,z,1);
        L2.push_back(pact); 
      }	  
  std::cout << "before insert !\n"; std::cout.flush();  	  
  T.insert(L2.begin(), L2.end());
  std::cout << "after insert !\n"; std::cout.flush();
  
  // check the result ...
  bool valid = T.is_valid(true);
  
  if (valid) std::cout << "valid result !\n";
  else std::cout << "NON - valid result !\n";
  

  return 0;
}
