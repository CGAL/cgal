#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Random.h>
#include <CGAL/point_generators_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

int main(int argc, char* argv[])
{
  int no_of_pts;
  double radius;
  if(argc > 1)
  {
    no_of_pts = pow(10, atoi(argv[1]));
    radius = atoi(argv[2]);
  }
  else
  {
    no_of_pts = 1000000;
    radius = 1;
  }

  CGAL::Random random(7);
  typedef CGAL::Creator_uniform_3<double, K::Point_3> Creator;

  std::cout << no_of_pts << std::endl;

  CGAL::Random_points_on_sphere_3<K::Point_3, Creator> on_sphere(radius, random);
  for(int count=0; count<no_of_pts; ++count)
  {
    K::Point_3 p = *on_sphere;
    ++on_sphere;

    std::cout << p << std::endl;
  }

  return EXIT_SUCCESS;
}
