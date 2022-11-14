#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/convex_hull_3.h>

#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
typedef K::Point_3                                               Point_3;
typedef CGAL::Surface_mesh<Point_3>                              Mesh;

int main(int argc, char* argv[])
{
  const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("meshes/star.off");
  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(filename, sm))
  {
    std::cerr<< "Cannot open input file." <<std::endl;
    return 1;
  }

  //This will contain the extreme vertices
  std::vector<Mesh::Vertex_index> extreme_vertices;

  //call the function with the traits adapter for vertices
  CGAL::extreme_points_3(vertices(sm), std::back_inserter(extreme_vertices),
                         CGAL::make_extreme_points_traits_adapter(sm.points()));

  //print the number of extreme vertices
  std::cout << "There are  " << extreme_vertices.size() << " extreme vertices in this mesh." << std::endl;

  return 0;
}
