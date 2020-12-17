#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_reconstruction_3.h>
#include <CGAL/IO/OFF_reader.h>
#include <CGAL/IO/PLY_writer.h>
#include <CGAL/Surface_mesh.h>

using SCF   = CGAL::Simple_cartesian<float>;
using SCD   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

using Kernel    = EPICK;
using Point_3   = typename Kernel::Point_3;
using Segment_3 = typename Kernel::Segment_3;

using Surface_mesh = CGAL::Surface_mesh<Point_3>;
using KSR = CGAL::Kinetic_shape_reconstruction_3<Kernel>;

struct Polygon_map {

  using key_type   = std::vector<std::size_t>;
  using value_type = std::vector<Point_3>;
  using reference  = value_type;
  using category   = boost::readable_property_map_tag;

  const std::vector<Point_3>& points;
  Polygon_map(
    const std::vector<Point_3>& vertices) :
  points(vertices)
  { }

  friend reference get(const Polygon_map& map, const key_type& face) {
    reference polygon;
    polygon.reserve(face.size());
    std::transform(
      face.begin(), face.end(),
      std::back_inserter(polygon),
      [&](const std::size_t vertex_index) -> Point_3 {
        return map.points[vertex_index];
      });
    return polygon;
  }
};

int main(const int argc, const char** argv) {

  // Input.
  std::string input_filename = (argc > 1 ? argv[1] : "data/test_1_polygon_a.off");
  std::ifstream input_file(input_filename);

  std::vector<Point_3> input_vertices;
  std::vector< std::vector<std::size_t> > input_faces;

  if (!CGAL::read_OFF(input_file, input_vertices, input_faces)) {
    std::cerr << "ERROR: can't read the file " << input_filename << "!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << std::endl;
  std::cout << "--- INPUT STATS: " << std::endl;
  std::cout << "* input kernel: " << boost::typeindex::type_id<Kernel>().pretty_name() << std::endl;
  std::cout << "* number of input vertices: " << input_vertices.size() << std::endl;
  std::cout << "* number of input faces: " << input_faces.size() << std::endl;

  // Algorithm.
  const bool debug   = true;
  const bool verbose = true;
  KSR ksr(verbose, debug);
  const unsigned int k = (argc > 2 ? std::atoi(argv[2]) : 1);
  std::cout << "* input k: " << k << std::endl;
  const unsigned int n = 0;
  const unsigned int num_blocks = std::pow(n + 1, 3);
  std::cout << "* input blocks: " << num_blocks << std::endl;
  const Polygon_map polygon_map(input_vertices);
  const bool is_success = ksr.partition(input_faces, polygon_map, k, n);
  assert(is_success);

  // Output.
  const int num_support_planes = ksr.number_of_support_planes();
  CGAL_assertion(num_support_planes > 6);

  // Vertices.
  const std::size_t num_vertices = ksr.number_of_vertices(-1);
  std::vector<Point_3> output_vertices;
  ksr.output_partition_vertices(
    std::back_inserter(output_vertices), -1);
  assert(num_vertices == output_vertices.size());

  // Edges.
  const std::size_t num_edges = ksr.number_of_edges(-1);
  std::vector<Segment_3> output_edges;
  ksr.output_partition_edges(
    std::back_inserter(output_edges), -1);
  assert(num_edges == output_edges.size());

  // Faces.
  const std::size_t num_faces = ksr.number_of_faces(-1);
  // output_vertices.clear();
  // std::vector< std::vector<std::size_t> > output_faces;
  // ksr.output_partition_faces(
  //   std::back_inserter(output_vertices),
  //   std::back_inserter(output_faces), -1);
  // assert(num_faces == output_faces.size());

  // const int num_volume_levels = ksr.number_of_volume_levels();
  // CGAL_assertion(num_volume_levels > 0);

  // Volumes.
  const std::size_t num_volumes = ksr.number_of_volumes(-1);
  // std::vector<Surface_mesh> output_volumes;
  // ksr.output_partition_volumes(
  //   std::back_inserter(output_volumes), -1);
  // assert(num_volumes == output_volumes.size());

  std::cout << std::endl;
  std::cout << "--- OUTPUT STATS: " << std::endl;
  std::cout << "* number of output vertices: " << num_vertices << std::endl;
  std::cout << "* number of output edges: "    << num_edges    << std::endl;
  std::cout << "* number of output faces: "    << num_faces    << std::endl;
  std::cout << "* number of output volumes: "  << num_volumes  << std::endl;

  // Export.
  std::cout << std::endl;
  std::cout << "--- EXPORT: " << std::endl;

  // Vertices.
  std::string output_filename = "partition-vertices.xyz";
  std::ofstream output_file_vertices(output_filename);
  output_file_vertices.precision(20);
  for (const auto& output_vertex : output_vertices)
    output_file_vertices << output_vertex << std::endl;
  output_file_vertices.close();
  std::cout << "* partition vertices exported successfully" << std::endl;

  // Edges.
  output_filename = "partition-edges.polylines.txt";
  std::ofstream output_file_edges(output_filename);
  output_file_edges.precision(20);
  for (const auto& output_edge : output_edges)
    output_file_edges << "2 " << output_edge << std::endl;
  output_file_edges.close();
  std::cout << "* partition edges exported successfully" << std::endl;

  // Faces.
  // output_filename = "partition-faces.ply";
  // std::ofstream output_file_faces(output_filename);
  // output_file_faces.precision(20);
  // if (!CGAL::write_PLY(output_file_faces, output_vertices, output_faces)) {
  //   std::cerr << "ERROR: can't write to the file " << output_filename << "!" << std::endl;
  //   return EXIT_FAILURE;
  // }
  // output_file_faces.close();
  // std::cout << "* partition faces exported successfully" << std::endl;

  // Volumes.
  // output_filename = "partition-volume-";
  // for (std::size_t i = 0; i < num_volumes; ++i) {
  //   const auto output_file = output_filename + std::to_string(i) + ".ply";
  //   std::ofstream output_file_volume(output_file);
  //   output_file_volume.precision(20);
  //   if (!CGAL::write_ply(output_file_volume, output_volumes[i])) {
  //     std::cerr << "ERROR: can't write to the file " << output_file << "!" << std::endl;
  //     return EXIT_FAILURE;
  //   }
  //   output_file_volume.close();
  // }
  // std::cout << "* partition volumes exported successfully" << std::endl;

  std::cout << std::endl << "3D KINETIC DONE!" << std::endl << std::endl;
  return EXIT_SUCCESS;
}
