#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_partition_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/IO/PLY.h>

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
using KSP          = CGAL::Kinetic_shape_partition_3<EPICK>;
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
  std::string input_filename = (argc > 1 ? argv[1] : "../data/test-4-rnd-polygons-4-6.off");
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
  const unsigned int k = (argc > 2 ? std::atoi(argv[2]) : 1); // intersections

  // Algorithm.
  KSP ksp(CGAL::parameters::verbose(true).debug(true));

  ksp.insert(input_vertices, input_faces);

  Timer timer;
  timer.start();
  ksp.initialize(CGAL::parameters::bbox_dilation_ratio(1.1).reorient_bbox(false));

  ksp.partition(k);

  timer.stop();
  const FT time = static_cast<FT>(timer.time());

  // Output.
  CGAL::Linear_cell_complex_for_combinatorial_map<3, 3> lcc;
  ksp.get_linear_cell_complex(lcc);

  std::vector<unsigned int> cells = { 0, 2, 3 }, count;
  count = lcc.count_cells(cells);

  std::cout << "For k = " << k << ":" << std::endl << " vertices: " << count[0] << std::endl << " faces: " << count[2] << std::endl << " volumes: " << count[3] << std::endl;

  std::cout << std::endl << "3D kinetic partition created in " << time << " seconds!" << std::endl << std::endl;
  return EXIT_SUCCESS;
}
