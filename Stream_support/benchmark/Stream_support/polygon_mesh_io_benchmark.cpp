#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>

#include <array>
#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef K::FT                                                   FT;
typedef K::Point_3                                              Point_3;
typedef CGAL::Surface_mesh<Point_3>                             Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  std::vector<Point_3> points;
  std::vector<std::array<std::size_t,3>> polygons;

  if (argc!=2)
  {
    std::cerr << "Usage: " << argv[0] << " input.XXX \n";
    std::cerr << "See https://doc.cgal.org/latest/Stream_support/index.html#title11 for the list of supported file formats.\n";
    return 0;
  }

  CGAL::Timer timer;
  timer.start();
  if (!CGAL::IO::read_polygon_soup(argv[1], points, polygons))
  {
    std::cerr << "Error reading " << argv[1] << "\n";
    return EXIT_FAILURE;
  }
  std::cout << "Read " << points.size() << " points and " << polygons.size() << " polygons in " << timer.time() << " seconds.\n";
  timer.stop();
  timer.reset();
  timer.start();

  Mesh mesh;
  PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);
  std::cout << "Constructed a polygon mesh with " << num_vertices(mesh) << " vertices and " << num_faces(mesh) << " faces in " << timer.time() << " seconds.\n";
  return 0;
}
