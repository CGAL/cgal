#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_reconstruction_3.h>
#include <CGAL/IO/OFF_reader.h>
#include <CGAL/IO/PLY_writer.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Real_timer.h>

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
using KSR          = CGAL::Kinetic_shape_reconstruction_3<Kernel>;
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

  // const FT x5 = 4.771745262374, y5 = 4.395963608911; // point J
  // const FT x6 = 9.000000000000, y6 = 5.000000000000; // point B
  // const FT xx = 7.423291494505, yy = 4.774755927786; // point Q
  // const FT x2 = 4.074202521868, y2 = 5.606021913821; // point M
  // const FT x3 = 13.82144413465, y3 = 3.186692216531; // point N

  // const Point_2 J(x5, y5);
  // const Point_2 B(x6, y6);
  // const Point_2 Q(xx, yy);
  // const Point_2 M(x2, y2);
  // const Point_2 N(x3, y3);

  // const FT a1 = x5-x6;
  // const FT b1 = y5-y6;
  // const FT c1 = x6*x6-x6*x5-y6*y5+y6*y6;

  // const FT d1 = (x6-x5)*(x6-x5)+(y6-y5)*(y6-y5);

  // const FT a2 = a1/d1;
  // const FT b2 = b1/d1;
  // const FT c2 = c1/d1;

  // const FT l1 = a2*xx+b2*yy+c2;
  // const FT l2 = FT(1)-l1;

  // const FT a3 = x2-x3;
  // const FT b3 = y2-y3;
  // const FT c3 = x3*x3-x3*x2-y3*y2+y3*y3;

  // const FT d2 = (x3-x2)*(x3-x2)+(y3-y2)*(y3-y2);

  // const FT a4 = a3/d2;
  // const FT b4 = b3/d2;
  // const FT c4 = c3/d2;

  // const FT m1 = a4*xx+b4*yy+c4;
  // const FT m2 = FT(1)-m1;

  // const FT a5 = x5*a2-x6*a2-x2*a4+x3*a4;
  // const FT b5 = x5*b2-x6*b2-x2*b4+x3*b4;
  // const FT c5 = x5*c2+x6-x6*c2-x2*c4-x3+x3*c4;

  // const FT a6 = y5*a2-y6*a2-y2*a4+y3*a4;
  // const FT b6 = y5*b2-y6*b2-y2*b4+y3*b4;
  // const FT c6 = y5*c2+y6-y6*c2-y2*c4-y3+y3*c4;

  // const FT x = (c5*b6-b5*c6)/(b5*a6-a5*b6);
  // const FT y = (-c5-a5*x)/(b5);

  // const FT lambda1 = a2*x+b2*y+c2;
  // const FT lambda2 = FT(1)-lambda1;

  // std::cout << "--debug--" << std::endl;
  // std::cout.precision(20);
  // std::cout << xx << " =? " << l1*x5+l2*x6 << std::endl;
  // std::cout << yy << " =? " << l1*y5+l2*y6 << std::endl;
  // std::cout << xx << " =? " << m1*x2+m2*x3 << std::endl;
  // std::cout << yy << " =? " << m1*y2+m2*y3 << std::endl;
  // std::cout << a5*xx+b5*yy+c5 << " =? " << 0 << std::endl;
  // std::cout << a6*xx+b6*yy+c6 << " =? " << 0 << std::endl;

  // std::cout << "--result--" << std::endl;
  // std::cout.precision(20);
  // std::cout << xx << " =? " << x << std::endl;
  // std::cout << yy << " =? " << y << std::endl;
  // std::cout << "lambda1 = " << lambda1 << std::endl;
  // std::cout << "lambda2 = " << lambda2 << std::endl;

  // exit(EXIT_SUCCESS);

  // Input.
  std::cout.precision(20);
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

  Timer timer;
  timer.start();
  const bool is_ksr_success = ksr.partition(
    input_faces, polygon_map, CGAL::parameters::
    k_intersections(k).
    n_subdivisions(subdiv).
    enlarge_bbox_ratio(eratio).
    reorient(orient));
  assert(is_ksr_success);
  const std::string success = is_ksr_success ? "SUCCESS" : "FAILED";
  timer.stop();
  const FT time = static_cast<FT>(timer.time());

  // Output.
  const int support_plane_idx = -1;
  const int num_support_planes = ksr.number_of_support_planes();
  assert(num_support_planes > 6);
  assert(ksr.support_plane_index(0) == 6);

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
  assert(num_volume_levels > 0);

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
    assert(sp_mesh.number_of_vertices() == ksr.number_of_vertices(i));
    assert(sp_mesh.number_of_edges()    == ksr.number_of_edges(i));
    assert(sp_mesh.number_of_faces()    == ksr.number_of_faces(i));
    support_planes.push_back(sp_mesh);
  }
  assert(support_planes.size() == num_support_planes);

  std::cout << std::endl;
  std::cout << "--- OUTPUT STATS: " << std::endl;
  std::cout << "* number of vertices: "       << num_vertices       << std::endl;
  std::cout << "* number of edges: "          << num_edges          << std::endl;
  std::cout << "* number of faces: "          << num_faces          << std::endl;
  std::cout << "* number of volumes: "        << num_volumes        << std::endl;
  std::cout << "* number of support planes: " << num_support_planes << std::endl;
  std::cout << "* number of volume levels: "  << num_volume_levels  << std::endl;

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

  std::cout << std::endl << "3D KINETIC " << success <<
  " in " << time << " seconds!" << std::endl << std::endl;
  return EXIT_SUCCESS;
}
