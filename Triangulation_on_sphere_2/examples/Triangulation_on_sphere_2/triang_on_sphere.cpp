#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_on_sphere_traits_2.h>
#include <CGAL/Delaunay_triangulation_on_sphere_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Delaunay_triangulation_on_sphere_traits_2<K>  Traits;
typedef CGAL::Delaunay_triangulation_on_sphere_2<Traits>    DToS2;

typedef Traits::Point_3                                     Point_3;

int main(int, char**)
{
  std::vector<Point_3> points;
  points.emplace_back( 2, 1, 1);
  points.emplace_back(-2, 1, 1); // not on the sphere
  points.emplace_back( 0, 1, 1);
  points.emplace_back( 1, 2, 1);
  points.emplace_back( 0, 1, 1); // duplicate of #3
  points.emplace_back( 1, 0, 1);
  points.emplace_back( 1, 1, 2);

  Traits traits(Point_3(1, 1, 1), 1); // sphere center on (1,1,1), with radius 1
  DToS2 dtos(traits);

  for(const Point_3& pt : points)
  {
    std::cout << "Inserting (" << pt
              << ") at squared distance " << CGAL::squared_distance(pt, traits.center())
              << " from the center of the sphere; is it on there sphere? "
              << (traits.is_on_sphere(pt) ? "yes" : "no") << std::endl;
    dtos.insert(pt);

    std::cout << "After insertion, the dimension of the triangulation is: " << dtos.dimension() << "\n";
    std::cout << "It has:\n";
    std::cout << dtos.number_of_vertices() << " vertices\n";
    std::cout << dtos.number_of_edges() << " edges\n";
    std::cout << dtos.number_of_solid_faces() << " solid faces\n";
    std::cout << dtos.number_of_ghost_faces() << " ghost faces\n" << std::endl;
  }

  CGAL::IO::write_OFF("result.off", dtos, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
