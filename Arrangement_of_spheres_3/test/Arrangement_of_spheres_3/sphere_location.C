#define CGAL_CHECK_EXACTNESS
#define CGAL_CHECK_EXPENSIVE

#include <CGAL/Arrangement_of_spheres_traits_3.h>



typedef Arrangement_of_spheres_traits_3 Traits;
typedef Traits::Sphere_point_3 SPoint_3;
typedef Traits::Event_point_3 EPoint_3;
typedef Traits::Point_3 Point;
typedef Traits::Vector_3 Vector;
typedef Traits::Sphere_3 Sphere;
typedef Traits::Line_3 Line;
typedef Traits::Plane_3 Plane;

void test(Traits::Sphere_location sl, Sphere s, int ans) {
  int l= sl(s);
  if (l != ans) {
    std::cerr << "Location failed for " << s << std::endl;
    std::cerr << "Got " << Traits::Sphere_location::decode(l)
	      << " expected " << Traits::Sphere_location::decode(ans) << std::endl;


  }
  
}

int main(int, char *[]){
  
  Traits tr;

  EPoint_3 ep(Sphere(Point(0,0,4), 16), Line(Point(0,0,0), Vector(0,0,1)));
  
  Traits::Sphere_location sl= tr.sphere_location_object(ep);
  // enum Location {L_BIT=1, R_BIT=2, T_BIT=4, B_BIT=8, IN_BIT=16 };
  std::cout << sl.sp_ << std::endl;
  test(sl, Sphere(Point(0,0,0), 1), 31 );
  test(sl, Sphere(Point(1,1,0), 1), 9);
  test(sl, Sphere(Point(1,0,0), 1), 29);
  test(sl, Sphere(Point(-1,0,0), 1), 30);
  test(sl, Sphere(Point(0,1,0), 1), 27);
  test(sl, Sphere(Point(0,-1,0), 1),23);
  test(sl, Sphere(Point(0,-1,1), 1), 7);

  return EXIT_SUCCESS;
}
