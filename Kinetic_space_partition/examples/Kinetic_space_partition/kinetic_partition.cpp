#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_space_partition_3.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/polygon_soup_io.h>

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

using Kernel     = EPICK;
using FT         = typename Kernel::FT;
using Point_3    = typename Kernel::Point_3;

using Surface_mesh = CGAL::Surface_mesh<Point_3>;
using KSP          = CGAL::Kinetic_space_partition_3<EPICK>;
using Timer        = CGAL::Real_timer;

int main(int argc, char** argv)
{
  // Reading polygons from file
  std::string input_filename = (argc > 1 ? argv[1] : "data/test-4-rnd-polygons-4-6.off");
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

  std::cout << "--- INPUT STATS: \n* number of polygons: " << input_faces.size() << std::endl;

  // Parameters.
  const unsigned int k = (argc > 2 ? std::atoi(argv[2]) : 1);

  // Initialization of Kinetic_space_partition_3 object.
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
  typedef CGAL::Linear_cell_complex_traits<3, EPECK> LCC_Traits;
  CGAL::Linear_cell_complex_for_combinatorial_map<3, 3, LCC_Traits, typename KSP::Linear_cell_complex_min_items> lcc;
  ksp.get_linear_cell_complex(lcc);

  std::vector<unsigned int> cells = { 0, 2, 3 }, count;
  count = lcc.count_cells(cells);

  std::cout << "For k = " << k << ":\n" << " vertices: " << count[0] << "\n faces: " << count[2] << "\n volumes: " << count[3] << std::endl;

  std::cout << "\n3D kinetic partition created in " << time << " seconds!" << std::endl;

  return EXIT_SUCCESS;
}
