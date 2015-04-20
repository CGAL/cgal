#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/spatial_sort_on_sphere.h>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef K::Point_3                                              Point;
typedef K::Vector_3                                             Vector;
typedef K::Sphere_3                                             Sphere;
typedef CGAL::Creator_uniform_3<double,Point>                   Creator_3;

int main ()
{
  std::size_t size = 32;
  CGAL::Random random (42);
  std::vector<Point> v;

  // unit sphere
  std::cout << "UNIT SPHERE: " << std::endl;

  v.reserve(size);
                                                                    
  CGAL::Random_points_on_sphere_3<Point> unit_sphere(1.0, random);  // generate points
  for (std::size_t i = 0; i < size; ++i) v.push_back(*unit_sphere++);

  CGAL::spatial_sort_on_sphere(v.begin(), v.end());                 // sort

  for(std::size_t i=0; i<size; ++i) std::cout << v[i] << std::endl; //output

  v.clear();

  // given sphere
  std::cout << "GIVEN SPHERE: " << std::endl;

  v.reserve(size);

  Vector trans = Vector(3,4,5);
  Sphere sphere = Sphere(CGAL::ORIGIN + trans, 4);
  CGAL::Random_points_on_sphere_3<Point> given_sphere(2.0, random);  // generate points
  for (std::size_t i = 0; i < size; ++i) v.push_back(*given_sphere++ + trans);

  CGAL::spatial_sort_on_sphere(v.begin(), v.end(),                   // sort
    sphere.squared_radius(), sphere.center());

  for(std::size_t i=0; i<size; ++i) std::cout << v[i] << std::endl; //output

  return 0;
}
