#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_partitioning_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/IO/PLY.h>

using SCF   = CGAL::Simple_cartesian<float>;
using SCD   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

using Traits = typename CGAL::Kinetic_shape_partitioning_traits_3<EPICK, EPECK, std::vector<typename EPICK::Point_3>, CGAL::Identity_property_map<typename EPICK::Point_3> >;

using Kernel     = EPICK;
using FT         = typename Kernel::FT;
using Point_2    = typename Kernel::Point_2;
using Point_3    = typename Kernel::Point_3;
using Segment_3  = typename Kernel::Segment_3;
using Triangle_2 = typename Kernel::Triangle_2;

using Surface_mesh = CGAL::Surface_mesh<Point_3>;
using KSR          = CGAL::Kinetic_shape_partitioning_3<Traits>;
using Timer        = CGAL::Real_timer;

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
  std::cout.precision(20);
  const auto kernel_name = boost::typeindex::type_id<Kernel>().pretty_name();
  std::string input_filename = (argc > 1 ? argv[1] : "data/stress-test-0/test-1-polygon-a.off");
  std::ifstream input_file_off(input_filename);
  std::ifstream input_file_ply(input_filename);

  std::vector<Point_3> input_vertices;
  std::vector< std::vector<std::size_t> > input_faces;

  if (CGAL::IO::read_OFF(input_file_off, input_vertices, input_faces)) {
    std::cout << "* reading the OFF file: " << input_filename << "!" << std::endl;
    input_file_off.close();
  } else if (CGAL::IO::read_PLY(input_file_ply, input_vertices, input_faces)) {
    std::cout << "* reading the PLY file: " << input_filename << "!" << std::endl;
    input_file_ply.close();
  } else {
    std::cerr << "ERROR: can't read the OFF/PLY file " << input_filename << "!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << std::endl;
  std::cout << "--- INPUT STATS: " << std::endl;
  std::cout << "* used kernel: "        << kernel_name        << std::endl;
  std::cout << "* number of polygons: " << input_faces.size() << std::endl;

  // Parameters.
  const bool         verbose = true;
  const bool         debug   = false;
  const unsigned int k = (argc > 2 ? std::atoi(argv[2]) : 1); // intersections
  const double       eratio  = 1.1;
  const bool         orient  = false;

  // Algorithm.
  KSR ksr(verbose, debug);
  const Polygon_map polygon_map(input_vertices);

  Timer timer;
  timer.start();
  bool is_ksr_success = ksr.initialize(
    input_faces, polygon_map, CGAL::parameters::
    bbox_dilation_ratio(eratio).
    reorient_bbox(orient));

  if (is_ksr_success)
    is_ksr_success = ksr.partition(k);

  assert(is_ksr_success);
  const std::string success = is_ksr_success ? "SUCCESS" : "FAILED";
  timer.stop();
  const FT time = static_cast<FT>(timer.time());

  // Output.
  const int support_plane_idx = -1;
  const int num_support_planes = ksr.number_of_support_planes();
  assert(num_support_planes > 6);
  assert(ksr.support_plane_index(0) == 6);/*

  // Vertices.
  const std::size_t num_vertices = ksr.number_of_vertices(support_plane_idx);
  std::vector<Point_3> output_vertices;
  ksr.output_partition_vertices(
    std::back_inserter(output_vertices), support_plane_idx);
  assert(num_vertices == output_vertices.size());

  // Edges.
  const std::size_t num_edges = ksr.number_of_edges(support_plane_idx);
  std::vector<Segment_3> output_edges;
  ksr.output_partition_edges(
    std::back_inserter(output_edges), support_plane_idx);
  assert(num_edges == output_edges.size());

  // Faces.
  const std::size_t num_faces = ksr.number_of_faces(support_plane_idx);
  std::vector< std::vector<std::size_t> > output_faces;
  ksr.output_partition_faces(
    std::back_inserter(output_faces), support_plane_idx, 6);
  assert(num_faces >= output_faces.size());

  // Volumes.
  const std::size_t num_volumes = ksr.number_of_volumes();
  std::vector<Surface_mesh> output_volumes;
  ksr.output_partition_volumes(
    std::back_inserter(output_volumes));
  assert(num_volumes == output_volumes.size());

  // Support planes.
  std::vector<Surface_mesh> support_planes;
  support_planes.reserve(num_support_planes);
  for (int i = 0; i < num_support_planes; ++i) {
    Surface_mesh sp_mesh;
    ksr.output_support_plane(sp_mesh, i);
    assert(sp_mesh.number_of_vertices() == ksr.number_of_vertices(i));
    assert(sp_mesh.number_of_edges()    == ksr.number_of_edges(i));
    assert(sp_mesh.number_of_faces()    == ksr.number_of_faces(i));
    support_planes.push_back(sp_mesh);
  }
  assert(support_planes.size() == static_cast<std::size_t>(num_support_planes));

  std::cout << std::endl;
  std::cout << "--- OUTPUT STATS: " << std::endl;
  std::cout << "* number of vertices: "       << num_vertices       << std::endl;
  std::cout << "* number of edges: "          << num_edges          << std::endl;
  std::cout << "* number of faces: "          << num_faces          << std::endl;
  std::cout << "* number of volumes: "        << num_volumes        << std::endl;
  std::cout << "* number of support planes: " << num_support_planes << std::endl;
  std::cout << "* number of events: "         << num_events         << std::endl;

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
  output_filename = "partition-faces.ply";
  std::ofstream output_file_faces(output_filename);
  output_file_faces.precision(20);
  if (!CGAL::IO::write_PLY(output_file_faces, output_vertices, output_faces)) {
    std::cerr << "ERROR: can't write to the file " << output_filename << "!" << std::endl;
    return EXIT_FAILURE;
  }
  output_file_faces.close();
  std::cout << "* partition faces exported successfully" << std::endl;

  // Volumes.
  output_filename = "partition-volume-";
  for (std::size_t i = 0; i < num_volumes; ++i) {
    const auto output_file = output_filename + std::to_string(i) + ".ply";
    std::ofstream output_file_volume(output_file);
    output_file_volume.precision(20);
    if (!CGAL::IO::write_PLY(output_file_volume, output_volumes[i])) {
      std::cerr << "ERROR: can't write to the file " << output_file << "!" << std::endl;
      return EXIT_FAILURE;
    }
    output_file_volume.close();
  }
  std::cout << "* partition volumes exported successfully" << std::endl;

  // Support planes.
  // For big data sets, this one will print too many data,
  // so we comment it out.
  // for (std::size_t i = 0; i < support_planes.size(); ++i) {
  //   const std::string filename = "support_plane-" + std::to_string(i) + ".ply";
  //   std::ofstream output_file_support_plane(filename);
  //   output_file_support_plane.precision(20);
  //   CGAL::IO::write_PLY(output_file_support_plane, support_planes[i]);
  //   output_file_support_plane.close();
  // }
  // std::cout << "* partition support planes exported successfully" << std::endl;*/

  std::cout << std::endl << "3D KINETIC " << success <<
  " in " << time << " seconds!" << std::endl << std::endl;
  return EXIT_SUCCESS;
}
