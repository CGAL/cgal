#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_reconstruction_3.h>
#include <CGAL/IO/OFF_reader.h>
#include <CGAL/IO/PLY_writer.h>

using SCF   = CGAL::Simple_cartesian<float>;
using SCD   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

using Kernel    = EPECK;
using Point_3   = typename Kernel::Point_3;
using Segment_3 = typename Kernel::Segment_3;

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

int main (int argc, char** argv) {

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
  KSR ksr;
  const unsigned int k = (argc > 2 ? std::atoi(argv[2]) : 1);
  std::cout << "* input k: " << k << std::endl;
  Polygon_map polygon_map(input_vertices);
  const bool is_success = ksr.partition(input_faces, polygon_map, k);
  assert(is_success);

  // Output.
  std::vector<Segment_3> output_edges;
  // ksr.output_partition_edges_to_segment_soup(std::back_inserter(output_edges));

  std::vector<Point_3> output_vertices;
  std::vector< std::vector<std::size_t> > output_faces;
  // ksr.output_partition_faces_to_polygon_soup(
  //   std::back_inserter(output_vertices), std::back_inserter(output_faces));

  std::cout << std::endl;
  std::cout << "--- OUTPUT STATS: " << std::endl;
  std::cout << "* number of output edges: " << output_edges.size() << std::endl;
  std::cout << "* number of output vertices: " << output_vertices.size() << std::endl;
  std::cout << "* number of output faces: " << output_faces.size() << std::endl;

  // Export.
  std::cout << std::endl;
  std::cout << "--- EXPORT: " << std::endl;

  std::string output_filename = "partition-edges.polylines";
  std::ofstream output_file_edges(output_filename);
  output_file_edges.precision(20);
  for (const auto& output_edge : output_edges)
    output_file_edges << "2 " << output_edge << std::endl;
  output_file_edges.close();
  std::cout << "* edges exported successfully" << std::endl;

  output_filename = "partition-faces.ply";
  std::ofstream output_file_faces(output_filename);
  output_file_faces.precision(20);
  if (!CGAL::write_PLY(output_file_faces, output_vertices, output_faces)) {
    std::cerr << "ERROR: can't write to the file " << output_filename << "!" << std::endl;
    return EXIT_FAILURE;
  }
  output_file_faces.close();
  std::cout << "* faces exported successfully" << std::endl;
  std::cout << std::endl << "3D KINETIC DONE!" << std::endl << std::endl;

  return EXIT_SUCCESS;
}
