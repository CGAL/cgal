#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Convex_hull_3/predicates.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
typedef K::Point_3                                               Point_3;
typedef CGAL::Surface_mesh<Point_3>                              Mesh;

int main(int argc, char* argv[])
{
  const std::string f1 = (argc>1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");
  const std::string f2 = (argc>2) ? argv[2] : CGAL::data_file_path("meshes/sphere.off");
  Mesh sm1, sm2;
  if(!CGAL::IO::read_polygon_mesh(f1, sm1))
  {
    std::cerr<< "Cannot open " << f1 <<std::endl;
    return 1;
  }
  if(!CGAL::IO::read_polygon_mesh(f2, sm2))
  {
    std::cerr<< "Cannot open " << f2 <<std::endl;
    return 1;
  }

  //This will contain the extreme vertices
  std::vector<Mesh::Vertex_index> extreme_vertices1, extreme_vertices2;

  //call the function with the traits adapter for vertices
  CGAL::extreme_points_3(vertices(sm1), std::back_inserter(extreme_vertices1),
                         CGAL::make_extreme_points_traits_adapter(sm1.points()));
  CGAL::extreme_points_3(vertices(sm2), std::back_inserter(extreme_vertices2),
                         CGAL::make_extreme_points_traits_adapter(sm2.points()));

  //print the number of extreme vertices
  std::cout << "There are " << extreme_vertices1.size() << " and "
                            << extreme_vertices2.size()
                            << " extreme vertices in the input meshes.\n";

  bool res = CGAL::Convex_hull_3::do_intersect(extreme_vertices1, extreme_vertices2,
                                               CGAL::parameters::point_map(sm1.points()),
                                               CGAL::parameters::point_map(sm2.points()));

  std::cout << "do convex hulls intersect? " << res << "\n";

  return 0;
}
