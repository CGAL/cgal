#include <vector>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/poisson_eliminate.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/IO/write_points.h>

using Point_3 = CGAL::Simple_cartesian<double>::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;

int main(int argc, char* argv[])
{
  std::vector<Point_3> points;
  std::size_t sampling_size = 20000;
  std::size_t target_size = 10000;
  Mesh mesh;

  if (argc == 4) {
    target_size = std::atoi(argv[2]);
    sampling_size = std::atoi(argv[3]);
    if (sampling_size < target_size) {
      std::cout << "usage: poisson_eliminate_from_mesh_example input_mesh target_size sampling_size"
                << std::endl << "with sampling_size > target_size, e.g., by a factor of 2" << std::endl;
      return -1;
    }
  }

  if (argc == 3) {
    target_size = std::atoi(argv[2]);
    sampling_size = 2 * target_size;
  }

  if (argc < 2)
    CGAL::IO::read_polygon_mesh(CGAL::data_file_path("meshes/elephant.off"), mesh);
  else
    CGAL::IO::read_polygon_mesh(argv[1], mesh);

  if (mesh.number_of_faces() == 0) {
    std::cout << "mesh could not be loaded" << std::endl;
    std::cout << "usage: poisson_eliminate_from_mesh_example input_mesh target_size sampling_size"
              << std::endl << "with sampling_size > target_size, e.g., by a factor of 2" << std::endl;
    return -1;
  }

  points.reserve(sampling_size);
  CGAL::Polygon_mesh_processing::sample_triangle_mesh(mesh,
    std::back_inserter(points),
    CGAL::parameters::number_of_points_on_faces(sampling_size)
    .do_sample_vertices(false)
    .do_sample_edges(false));

  CGAL::IO::write_points("out_sampled.xyz", points, CGAL::parameters::stream_precision(17));

  std::vector<Point_3> output;
  output.reserve(target_size);

  CGAL::poisson_eliminate(points, target_size, std::back_inserter(output));

  CGAL::IO::write_points("out_poisson_eliminated.xyz", output, CGAL::parameters::stream_precision(17));

  return 0;
}
