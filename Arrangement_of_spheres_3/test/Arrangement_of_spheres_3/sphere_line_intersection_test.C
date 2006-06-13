#include <CGAL/Arrangement_of_spheres_traits_3.h>
#include <CGAL/Arrangement_of_spheres_3/Filtered_sphere_line_intersection.h>

typedef Arrangement_of_spheres_traits_3 Traits;
typedef Traits::Sphere_point_3 SLI;
typedef Traits::Point_3 Point;
typedef Traits::Vector_3 Vector;
typedef Traits::Sphere_3 Sphere;
typedef Traits::Line_3 Line;
typedef Traits::Plane_3 Plane;


template <class T>
void test(){
  std::cout << "Test 1: " << std::endl;
  {
    T s0(Sphere(Point(0,0,0), 1), Line(Point(0, 1, 1), Vector(0,0,1)));
    T s1(Sphere(Point(0,0,0), 1.00001), Line(Point(0, 1, 1), Vector(0,0,1)));
    T s2(Sphere(Point(0,0,0), 1.00001), Line(Point(0, 1, 1), Vector(0,0,-1)));

    std::cout << s0 << std::endl;
    std::cout << s1 << std::endl;
    std::cout << s2 << std::endl;
    
    //CGAL::Oriented_side os= oriented_side(
    CGAL::Comparison_result c0= s0.compare(s1, typename T::Coordinate_index(0));
    CGAL::Comparison_result c1= s0.compare(s1, typename T::Coordinate_index(1));
    CGAL::Comparison_result c2= s0.compare(s1, typename T::Coordinate_index(2));
    CGAL::Comparison_result c3= s0.compare(s2, typename T::Coordinate_index(2));
    CGAL::Comparison_result c4= s1.compare(s2, typename T::Coordinate_index(2));
    CGAL_assertion(c0== CGAL::EQUAL);
    CGAL_assertion(c1== CGAL::EQUAL); 
    CGAL_assertion(c2== CGAL::LARGER);
    CGAL_assertion(c3== CGAL::SMALLER);
    CGAL_assertion(c4== CGAL::SMALLER);
  }
  std::cout << "Test 2: " << std::endl;
  {
    T s0(Sphere(Point(1,2,3), 1), Line(Point(1, 1.5, 2.5), Vector(1,-.5,-1)));
    T s1(Sphere(Point(1,2,3), 1.00001), Line(Point(1, 1.5, 2.5), Vector(-1,-.5,-.25)));
    T s2(Sphere(Point(1,2,3), 1.00001), Line(Point(1, 1.5, 2.5), Vector(1,.5,.25)));
    CGAL_assertion(s0.is_valid());
    CGAL_assertion(s1.is_valid());
    CGAL_assertion(s2.is_valid());
    std::cout << s0 << std::endl;
    std::cout << s1 << std::endl;
    std::cout << s2 << std::endl;
    
    //CGAL::Oriented_side os= oriented_side(
    CGAL::Comparison_result c0= s0.compare(s1, typename T::Coordinate_index(0));
    CGAL::Comparison_result c1= s0.compare(s1, typename T::Coordinate_index(1));
    CGAL::Comparison_result c2= s0.compare(s1, typename T::Coordinate_index(2));
    CGAL::Comparison_result c3= s0.compare(s2, typename T::Coordinate_index(2));
    CGAL::Comparison_result c4= s1.compare(s2, typename T::Coordinate_index(2));
    CGAL_assertion(c0== CGAL::SMALLER);
    CGAL_assertion(c1== CGAL::SMALLER); 
    CGAL_assertion(c2== CGAL::LARGER);
    CGAL_assertion(c3== CGAL::LARGER);
    CGAL_assertion(c4== CGAL::LARGER);
  }
  /*std::cout << "Test 3 " <<std::endl;
  {
    T s0(Sphere(Point(0,0,0), 1), Line(Point(0, 0, 0), Vector(1,1,1)));
    Plane p(Point(0,0,0), Vector(1,1,1));
    CGAL::Oriented_side os= oriented_side(p, s0);
    CGAL_assertion(os== CGAL::ON_NEGATIVE_SIDE);
    
    Plane op=p.opposite();
    CGAL::Oriented_side oos= oriented_side(op, s0);
    CGAL_assertion(oos== CGAL::ON_POSITIVE_SIDE);
    }*/
}

int main(int, char *[]){
  
  test<SLI>();
  test<Filtered_sphere_line_intersection<Traits::Geometric_traits, 2> >();
  
  return EXIT_SUCCESS;
}
