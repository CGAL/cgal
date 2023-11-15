#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_partition_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/polygon_soup_io.h>

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

double add_polys = 0, intersections = 0, iedges = 0, ifaces = 0, mapping = 0;

void get_affine_transformation(std::vector<Point_3>& pts, std::vector<std::vector<std::size_t> >& polys) {
  // Projection in 2d

  FT minz = (std::numeric_limits<FT>::max)(), maxz = -(std::numeric_limits<FT>::max)();
  std::vector<bool> used_pts(pts.size(), false);
  std::size_t count = 0;
  for (std::size_t i = 0; i < polys.size(); i++)
    for (std::size_t j = 0; j < polys[i].size(); j++) {
      if (!used_pts[polys[i][j]]) {
        used_pts[polys[i][j]] = true;
        count++;
      }
    }

  std::vector<Point_2> pts2d;
  pts2d.reserve(count);

  for (std::size_t i = 0; i < used_pts.size(); i++) {
    if (used_pts[i])
      pts2d.push_back(Point_2(pts[i].x(), pts[i].y()));
  }

  std::vector<Point_2> ch;
  CGAL::convex_hull_2(pts2d.begin(), pts2d.end(), std::back_inserter(ch));

  std::vector<Point_2> bbox;
  bbox.reserve(8);
  CGAL::min_rectangle_2(ch.begin(), ch.end(), std::back_inserter(bbox));
}

int main(const int argc, const char** argv) {

  // Reading polygons from file
  const auto kernel_name = boost::typeindex::type_id<Kernel>().pretty_name();
  std::string input_filename = (argc > 1 ? argv[1] : "../data/test-4-rnd-polygons-4-6.off");
  std::ifstream input_file(input_filename);

  std::vector<Point_3> input_vertices;
  std::vector<std::vector<std::size_t> > input_faces;

  if (CGAL::IO::read_polygon_soup(input_filename, input_vertices, input_faces)) {
    std::cout << "* reading the file: " << input_filename << "!" << std::endl;
    input_file.close();
  } else {
    std::cerr << "ERROR: can't read the file " << input_filename << "!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << std::endl;
  std::cout << "--- INPUT STATS: " << std::endl;
  std::cout << "* used kernel: "        << kernel_name        << std::endl;
  std::cout << "* number of polygons: " << input_faces.size() << std::endl;

  // Parameters.
  const unsigned int k = (argc > 2 ? std::atoi(argv[2]) : 1);

  // Initialization of Kinetic_shape_partition_3 object.
  // 'debug' set to true exports intermediate results into files in the working directory.
  // The resulting volumes are exported into a volumes folder, if the folder already exists.
  KSP ksp(CGAL::parameters::verbose(true).debug(true));

  // Providing input polygons.
  ksp.insert(input_vertices, input_faces);

  Timer timer;
  timer.start();

  // 'initialize' creates the intersection graph that is used for the partition.
  ksp.initialize(CGAL::parameters::bbox_dilation_ratio(1.1).reorient_bbox(false));

  // Creating the partition with allowing up to 'k' intersections for each kinetic polygon.
  ksp.partition(k);

  timer.stop();
  const FT time = static_cast<FT>(timer.time());

  // Access the kinetic partition via linear cell complex.
  typedef CGAL::Linear_cell_complex_traits<3, CGAL::Exact_predicates_exact_constructions_kernel> LCC_Traits;
  CGAL::Linear_cell_complex_for_combinatorial_map<3, 3, LCC_Traits, typename KSP::Linear_cell_complex_min_items> lcc;
  ksp.get_linear_cell_complex(lcc);

  std::vector<unsigned int> cells = { 0, 2, 3 }, count;
  count = lcc.count_cells(cells);

  std::cout << "For k = " << k << ":" << std::endl << " vertices: " << count[0] << std::endl << " faces: " << count[2] << std::endl << " volumes: " << count[3] << std::endl;

  std::cout << std::endl << "3D kinetic partition created in " << time << " seconds!" << std::endl << std::endl;
  return EXIT_SUCCESS;
}
