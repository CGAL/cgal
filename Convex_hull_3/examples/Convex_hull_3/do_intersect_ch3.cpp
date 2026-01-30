#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Convex_hull_3/distance.h>
#include <CGAL/Convex_hull_3/do_intersect.h>
#include <CGAL/Convex_hull_hierarchy.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
typedef K::Point_3                                               Point_3;
typedef CGAL::Surface_mesh<Point_3>                              Mesh;
typedef Mesh::Property_map<Mesh::Vertex_index, Point_3>          PointMap;
typedef CGAL::Convex_hull_hierarchy<Mesh>                        Convex_hull_hierarchy;

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

  //call the function with the traits adapter for vertices
  bool res = CGAL::Convex_hull_3::do_intersect(vertices(sm1), vertices(sm2),
                                               CGAL::parameters::point_map(sm1.points()),
                                               CGAL::parameters::point_map(sm2.points()));

  // Convex_hull_hierarchy is data structure for fast intersection test of Convex_hull
  Convex_hull_hierarchy hsm1(sm1);
  Convex_hull_hierarchy hsm2(sm2);
  res = CGAL::Convex_hull_3::do_intersect(hsm1, hsm2);
  std::cout << "do convex hulls intersect? " << std::boolalpha << res << "\n";


  return 0;
}