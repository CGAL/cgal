#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_on_sphere_2.h>
#include <CGAL/Projection_on_sphere_traits_3.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel          K;

typedef CGAL::Projection_on_sphere_traits_3<K>                       Traits;
typedef CGAL::Delaunay_triangulation_on_sphere_2<Traits>             DToS2;

typedef Traits::Point_3                                              Point_3;

int main(int, char**)
{
  std::vector<Point_3> points;
  points.emplace_back( 3,  1,  1);
  points.emplace_back(-8,  1,  1);
  points.emplace_back( 1,  2,  1);
  points.emplace_back( 1, -2,  1);
  points.emplace_back( 1,  1, 10);

  Traits traits(Point_3(4,1,1)); // radius is 1 by default
  DToS2 dtos(traits);

  Traits::Construct_point_on_sphere_2 cst = traits.construct_point_on_sphere_2_object();

  for(const auto& pt : points)
  {
    std::cout << "----- Inserting (" << pt
              << ") at squared distance " << CGAL::squared_distance(pt, traits.center())
              << " from the center of the sphere" << std::endl;
    dtos.insert(cst(pt));

    std::cout << "The triangulation now has dimension: " << dtos.dimension() << " and\n";
    std::cout << dtos.number_of_vertices() << " vertices" << std::endl;
    std::cout << dtos.number_of_edges() << " edges" << std::endl;
    std::cout << dtos.number_of_solid_faces() << " solid faces" << std::endl;
    std::cout << dtos.number_of_ghost_faces() << " ghost faces" << std::endl;
  }
  assert(dtos.is_valid());
  std::ofstream out("dtos.txt");
  out.precision(17);
  out << dtos << std::endl;
  out.close();

  CGAL::IO::write_OFF("dtos.off", dtos, CGAL::parameters::stream_precision(17));

  DToS2 dtos2(traits);
  std::ifstream in("dtos.txt");
  in >> dtos2;
  in.close();

    std::cout << "After write/read triangulation has dimension: " << dtos2.dimension() << " and\n";
    std::cout << dtos2.number_of_vertices() << " vertices" << std::endl;
    std::cout << dtos2.number_of_edges() << " edges" << std::endl;
    std::cout << dtos2.number_of_solid_faces() << " solid faces" << std::endl;
    std::cout << dtos2.number_of_ghost_faces() << " ghost faces" << std::endl;

  CGAL::IO::write_OFF("dtos2.off", dtos2, CGAL::parameters::stream_precision(17));


  return EXIT_SUCCESS;
}
