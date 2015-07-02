#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

#include <CGAL/point_generators_3.h>

#include <CGAL/Side_of_triangle_mesh.h>

#include <vector>
#include <fstream>
#include <limits>
#include <boost/foreach.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef OpenMesh::PolyMesh_ArrayKernelT< > Mesh;
typedef K::Point_3 Point;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

double max_coordinate(const Mesh& mesh)
{
  typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type VPmap;
  VPmap vpmap = get(CGAL::vertex_point,mesh);

  double max_coord = std::numeric_limits<double>::min();
  BOOST_FOREACH(vertex_descriptor v, vertices(mesh))
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
  const char* filename = (argc > 1) ? argv[1] : "data/eight.off";

  Mesh mesh;
  OpenMesh::IO::read_mesh(mesh, filename);
 
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

  std::cerr << "Total query size: " << points.size() << std::endl;
  std::cerr << "  " << nb_inside << " points inside " << std::endl;
  std::cerr << "  " << nb_boundary << " points on boundary " << std::endl;
  std::cerr << "  " << points.size() - nb_inside - nb_boundary << " points outside " << std::endl;

  return 0;
}
