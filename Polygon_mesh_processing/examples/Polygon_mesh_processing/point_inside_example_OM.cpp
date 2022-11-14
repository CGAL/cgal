#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <iostream>
#include <limits>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;

typedef OpenMesh::PolyMesh_ArrayKernelT< >                    Mesh;
typedef K::Point_3                                            Point;

typedef boost::graph_traits<Mesh>::vertex_descriptor          vertex_descriptor;

double max_coordinate(const Mesh& mesh)
{
  typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type VPmap;
  VPmap vpmap = get(CGAL::vertex_point,mesh);

  double max_coord = -std::numeric_limits<double>::infinity();
  for(vertex_descriptor v : vertices(mesh))
  {
    Point p = get(vpmap, v);
    max_coord = (std::max)(max_coord, p.x());
    max_coord = (std::max)(max_coord, p.y());
    max_coord = (std::max)(max_coord, p.z());
  }
  return max_coord;
}

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/eight.off");

  Mesh mesh;
  OpenMesh::IO::read_mesh(mesh, filename);
  if (CGAL::is_empty(mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::Side_of_triangle_mesh<Mesh, K> inside(mesh);

  double size = max_coordinate(mesh);

  unsigned int nb_points = 100;
  std::vector<Point> points;
  points.reserve(nb_points);
  CGAL::Random_points_in_cube_3<Point> gen(size);
  for (unsigned int i = 0; i < nb_points; ++i)
    points.push_back(*gen++);

  std::cout << "Test " << nb_points << " random points in cube "
    << "[-" << size << "; " << size <<"]" << std::endl;

  int nb_inside = 0;
  int nb_boundary = 0;
  for (std::size_t i = 0; i < nb_points; ++i)
  {
    CGAL::Bounded_side res = inside(points[i]);

    if (res == CGAL::ON_BOUNDED_SIDE) { ++nb_inside; }
    if (res == CGAL::ON_BOUNDARY) { ++nb_boundary; }
  }

  std::cout << "Total query size: " << points.size() << std::endl;
  std::cout << "  " << nb_inside << " points inside " << std::endl;
  std::cout << "  " << nb_boundary << " points on boundary " << std::endl;
  std::cout << "  " << points.size() - nb_inside - nb_boundary << " points outside " << std::endl;

  return EXIT_SUCCESS;
}
