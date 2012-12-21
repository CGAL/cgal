// example with random points on a sphere

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Origin.h>

#include <CGAL/surface_neighbor_coordinates_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                      Coord_type;
typedef K::Point_3                 Point_3;
typedef K::Vector_3                Vector_3;
typedef std::vector< std::pair< Point_3, K::FT  > >
                                   Point_coordinate_vector;

int main()
{

  int n=100;
  std::vector< Point_3> points;
  points.reserve(n);

  std::cout << "Generate " << n << " random points on a sphere."
	    << std::endl;
  CGAL::Random_points_on_sphere_3<Point_3> g(1);
  CGAL::cpp11::copy_n( g, n, std::back_inserter(points));

  Point_3 p(1, 0,0);
  Vector_3 normal(p-CGAL::ORIGIN);
  std::cout << "Compute surface neighbor coordinates for "
	    << p << std::endl;
  Point_coordinate_vector coords;
  CGAL::Triple< std::back_insert_iterator<Point_coordinate_vector>,
    K::FT, bool> result =
    CGAL::surface_neighbor_coordinates_3(points.begin(), points.end(),
					 p, normal,
					 std::back_inserter(coords),
					 K());
  if(!result.third){
    //Undersampling:
    std::cout << "The coordinate computation was not successful."
	      << std::endl;
    return 0;
  }
  K::FT norm = result.second;

  std::cout << "Testing the barycentric property " << std::endl;
  Point_3 b(0, 0,0);
  for(std::vector< std::pair< Point_3, Coord_type  > >::const_iterator
	it = coords.begin(); it!=coords.end(); ++it)
    b = b + (it->second/norm)* (it->first - CGAL::ORIGIN);

  std::cout <<"    weighted barycenter: " << b <<std::endl;
  std::cout << "    squared distance: " <<
    CGAL::squared_distance(p,b) <<std::endl;

  std::cout << "done" << std::endl;
  return 0;
}
