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

using Kernel     = EPICK;
using FT         = typename Kernel::FT;
using Point_2    = typename Kernel::Point_2;
using Point_3    = typename Kernel::Point_3;
using Segment_3  = typename Kernel::Segment_3;
using Triangle_2 = typename Kernel::Triangle_2;

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

  const FT x1 = 4.771745262374, y1 = 4.395963608911; // point J
  const FT x2 = 13.82144413465, y2 = 3.186692216531; // point N
  const FT xx = 7.423291494505, yy = 4.774755927786; // point Q
  const FT x3 = 6.166717106737, y3 = 5.086645983861; // point V
  const FT x4 = 4.074202521868, y4 = 5.606021913821; // point M
  const FT x5 = 2.000000000000, y5 = 4.000000000000; // point A
  const FT x6 = 9.000000000000, y6 = 5.000000000000; // point B

  const Point_2 J(x1, y1);
  const Point_2 N(x2, y2);
  const Point_2 Q(xx, yy);
  const Point_2 V(x3, y3);
  const Point_2 M(x4, y4);
  const Point_2 A(x5, y5);
  const Point_2 B(x6, y6);

  std::cout.precision(20);

  // std::cout << xx << " =? " << x << std::endl;
  // std::cout << yy << " =? " << y << std::endl;

  exit(EXIT_SUCCESS);

  // Input.
  const auto kernel_name = boost::typeindex::type_id<Kernel>().pretty_name();
  std::string input_filename = (argc > 1 ? argv[1] : "data/stress-test-0/test-1-polygon-a.off");
  std::ifstream input_file(input_filename);

  std::vector<Point_3> input_vertices;
  std::vector< std::vector<std::size_t> > input_faces;

  if (!CGAL::read_OFF(input_file, input_vertices, input_faces)) {
    std::cerr << "ERROR: can't read the file " << input_filename << "!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << std::endl;
  std::cout << "--- INPUT STATS: " << std::endl;
  std::cout << "* used kernel: "        << kernel_name        << std::endl;
  std::cout << "* number of polygons: " << input_faces.size() << std::endl;

  // Parameters.
  const bool         verbose = true;
  const bool         debug   = true;
  const unsigned int k = (argc > 2 ? std::atoi(argv[2]) : 1); // intersections
  const unsigned int subdiv  = 0;
  const double       eratio  = 1.1;
  const bool         orient  = false;

  // Algorithm.
  KSR ksr(verbose, debug);
  const Polygon_map polygon_map(input_vertices);
  const bool is_success = ksr.partition(
    input_faces, polygon_map, CGAL::parameters::
    k_intersections(k).
    n_subdivisions(subdiv).
    enlarge_bbox_ratio(eratio).
    reorient(orient));
  assert(is_success);

  // Output.
  const int support_plane_idx = -1;
  const int num_support_planes = ksr.number_of_support_planes();
  CGAL_assertion(num_support_planes > 6);
  CGAL_assertion(ksr.support_plane_index(0) == 6);

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
    std::back_inserter(output_faces), support_plane_idx);
  assert(num_faces == output_faces.size());

  int volume_level = -1;
  const int num_volume_levels = ksr.number_of_volume_levels();
  CGAL_assertion(num_volume_levels > 0);

  // Volumes.
  const std::size_t num_volumes = ksr.number_of_volumes(volume_level);
  std::vector<Surface_mesh> output_volumes;
  ksr.output_partition_volumes(
    std::back_inserter(output_volumes), volume_level);
  assert(num_volumes == output_volumes.size());

  // Support planes.
  std::vector<Surface_mesh> support_planes;
  support_planes.reserve(num_support_planes);
  for (int i = 0; i < num_support_planes; ++i) {
    Surface_mesh sp_mesh;
    ksr.output_support_plane(sp_mesh, i);
    CGAL_assertion(sp_mesh.number_of_vertices() == ksr.number_of_vertices(i));
    CGAL_assertion(sp_mesh.number_of_edges()    == ksr.number_of_edges(i));
    CGAL_assertion(sp_mesh.number_of_faces()    == ksr.number_of_faces(i));
    support_planes.push_back(sp_mesh);
  }
  CGAL_assertion(support_planes.size() == num_support_planes);

  std::cout << std::endl;
  std::cout << "--- OUTPUT STATS: " << std::endl;
  std::cout << "* number of vertices: "       << num_vertices       << std::endl;
  std::cout << "* number of edges: "          << num_edges          << std::endl;
  std::cout << "* number of faces: "          << num_faces          << std::endl;
  std::cout << "* number of volumes: "        << num_volumes        << std::endl;
  std::cout << "* number of support planes: " << num_support_planes << std::endl;

  // Export.
  std::cout << std::endl;
  std::cout << "--- EXPORT: " << std::endl;

  // Vertices.
  // std::string output_filename = "partition-vertices.xyz";
  // std::ofstream output_file_vertices(output_filename);
  // output_file_vertices.precision(20);
  // for (const auto& output_vertex : output_vertices)
  //   output_file_vertices << output_vertex << std::endl;
  // output_file_vertices.close();
  // std::cout << "* partition vertices exported successfully" << std::endl;

  // Edges.
  // output_filename = "partition-edges.polylines.txt";
  // std::ofstream output_file_edges(output_filename);
  // output_file_edges.precision(20);
  // for (const auto& output_edge : output_edges)
  //   output_file_edges << "2 " << output_edge << std::endl;
  // output_file_edges.close();
  // std::cout << "* partition edges exported successfully" << std::endl;

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

  // Support planes.
  // for (std::size_t i = 0; i < support_planes.size(); ++i) {
  //   const std::string filename = "support_plane-" + std::to_string(i) + ".ply";
  //   std::ofstream output_file_support_plane(filename);
  //   output_file_support_plane.precision(20);
  //   CGAL::write_ply(output_file_support_plane, support_planes[i]);
  //   output_file_support_plane.close();
  // }
  // std::cout << "* partition support planes exported successfully" << std::endl;

  std::cout << std::endl << "3D KINETIC DONE!" << std::endl << std::endl;
  return EXIT_SUCCESS;
}
