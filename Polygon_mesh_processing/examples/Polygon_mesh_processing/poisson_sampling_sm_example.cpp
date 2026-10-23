#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/distance.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <iostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef K::Point_3                                            Point;

typedef CGAL::Surface_mesh<Point>                             Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor            face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/eight.off");

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  const double sampling_radius = (argc > 2) ? std::atof(argv[2]) : 0.009;

  const std::size_t number_of_darts = (argc > 3) ? std::atof(argv[3]) : 20;

  CGAL::Random random_seed = (argc > 4) ? std::atof(argv[4]) : CGAL::get_default_random();

  std::vector<Point> points;
  PMP::sample_triangle_mesh(mesh,
                            std::back_inserter(points),
                            CGAL::parameters::use_poisson_disk_sampling_euclidean(true)
                                             .sampling_radius(sampling_radius)
                                             .number_of_darts(number_of_darts)
                                             .random_seed(random_seed));

  std::ofstream out("sampling.xyz");
  out << std::setprecision(17);
  std::copy(points.begin(), points.end(), std::ostream_iterator<Point>(out, "\n"));

  std::cout << points.size() << std::endl;
  return 0;
}
