#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/Kinetic_space_partition_3.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/IO/PLY.h>
#include <CGAL/Real_timer.h>

using SCD   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

std::size_t different = 0;
std::size_t no_input = 0;
std::size_t failed = 0;
std::size_t passed = 0;

template<typename Kernel, typename IntersectionKernel>
bool run_test(
  const std::string input_filename,
  const std::vector<unsigned int>& ks,
  const std::vector<std::vector<unsigned int> >& results) {

  using Point_3   = typename Kernel::Point_3;
  using KSP = CGAL::Kinetic_space_partition_3<Kernel, IntersectionKernel>;

  std::string filename = input_filename;
  std::ifstream input_file_off(filename);
  std::ifstream input_file_ply(filename);
  std::vector<Point_3> input_vertices;
  std::vector< std::vector<std::size_t> > input_faces;

  if (CGAL::IO::read_OFF(input_file_off, input_vertices, input_faces))
    input_file_off.close();
  else if (CGAL::IO::read_PLY(input_file_ply, input_vertices, input_faces))
    input_file_ply.close();
  else {
    std::cerr << "ERROR: can't read the OFF/PLY file " << filename << "!" << std::endl;
    no_input++;
    return false;
  }

  for (std::size_t i = 0; i < ks.size(); i++) {
    CGAL::Real_timer timer;
    timer.start();
    KSP ksp(CGAL::parameters::verbose(false).debug(false).bbox_dilation_ratio(1.004).max_octree_node_size(10));

    ksp.insert(input_vertices, input_faces, CGAL::parameters::verbose(false).debug(false).bbox_dilation_ratio(1.004));

    ksp.initialize(CGAL::parameters::verbose(false).debug(false).bbox_dilation_ratio(1.004).max_octree_node_size(10));
    //std::cout << std::endl << "--INPUT K: " << k << std::endl;
    ksp.partition(ks[i]);

    CGAL::Linear_cell_complex_for_combinatorial_map<3, 3, CGAL::Linear_cell_complex_traits<3, Kernel>, KSP::Linear_cell_complex_min_items> lcc;
    ksp.get_linear_cell_complex(lcc);

    //std::cout << "Time: " << timer.time() << " seconds." << std::endl;

    std::vector<unsigned int> cells = { 0, 2, 3 }, count;
    count = lcc.count_cells(cells);

    //std::cout << "v: " << count[0] << " f: " << count[2] << " v: " << count[3] << " for k: " << ks[i] << std::endl;

    if (ksp.number_of_volumes() == 0)
      failed++;
    else if (results[i][0] != count[0] || results[i][1] != count[2] || results[i][2] != count[3]) {
      std::cout << input_filename << std::endl;
      std::cout << "TEST differs: Partitioning has not expected number of vertices, faces or volumes for k = " << ks[i] << std::endl;

      std::cout << "Expectation:" << std::endl;
      std::cout << "v: " << results[i][0] << " f : " << results[i][1] << " v : " << results[i][2] << std::endl;
      std::cout << "Result k = " << " vertices : " << count[0] << " faces : " << count[2] << " volumes : " << count[3] << std::endl;
      std::cout << "v: " << count[0] << " f: " << count[2] << " v: " << count[3] << " for k: " << ks[i] << std::endl;

      different++;
    }
    else {
      passed++;
      //std::cout << "TEST PASSED k = " << ks[i] << " " << input_filename << std::endl;
    }
  }

  return true;
}

template<typename Kernel, typename IntersectionKernel>
void run_all_tests() {
  different = 0;
  std::cout.precision(10);
  std::vector< std::vector<double> > all_times;

  std::vector<std::vector<unsigned int> > results(4);

  CGAL::Real_timer timer;
  timer.start();

  results[0] = { 38, 46, 9 };
  results[1] = { 47, 67, 15 };
  run_test<Kernel, IntersectionKernel>("../data/edge-case-test/test-same-time.off", { 0, 1 }, results);

  // Edge tests.
  results[0] = { 18, 20, 4 };
  run_test<Kernel, IntersectionKernel>("../data/edge-case-test/test-2-polygons.off", { 1 }, results);

  results[0] = { 22, 25, 5 };
  run_test<Kernel, IntersectionKernel>("../data/edge-case-test/test-4-polygons.off", { 1 }, results);
  run_test<Kernel, IntersectionKernel>("../data/edge-case-test/test-5-polygons.off", { 1 }, results);

  results[0] = { 49, 71, 16 };
  run_test<Kernel, IntersectionKernel>("../data/edge-case-test/test-local-global-1.off", { 1 }, results);
  run_test<Kernel, IntersectionKernel>("../data/edge-case-test/test-local-global-2.off", { 1 }, results);

  // Stress tests 0.
  results[0] = { 14, 13, 2 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-1-polygon-a.off", { 0 }, results);
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-1-polygon-b.off", { 0 }, results);
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-1-polygon-c.off", { 0 }, results);
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-1-polygon-d.off", { 0 }, results);
  results[0] = { 19, 19, 3 };
  results[1] = { 20, 22, 4 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-2-polygons-ab.off", { 0, 1 }, results);
  results[0] = { 18, 18, 3 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-2-polygons-ac.off", { 0 }, results);
  results[0] = { 18, 18, 3 };
  results[1] = { 19, 21, 4 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-2-polygons-ad.off", { 0, 1 }, results);
  results[0] = { 18, 18, 3 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-2-polygons-bc.off", { 0 }, results);
  results[0] = { 18, 18, 3 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-2-polygons-bd.off", { 0 }, results);
  results[0] = { 19, 21, 4 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-2-polygons-cd.off", { 0 }, results);
  results[0] = { 24, 25, 4 };
  results[1] = { 25, 28, 5 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-3-polygons-abc.off", { 0, 1 }, results);
  results[0] = { 25, 26, 4 };
  results[1] = { 27, 32, 6 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-3-polygons-abd.off", { 0, 1 }, results);
  results[0] = { 24, 27, 5 };
  results[1] = { 26, 31, 6 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-3-polygons-acd.off", { 0, 1 }, results);
  results[0] = { 25, 28, 5 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-3-polygons-bcd.off", { 0 }, results);
  results[0] = { 31, 35, 6 };
  results[1] = { 34, 42, 8 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-4-polygons-abcd.off", { 0, 1 }, results);
  results[0] = { 47, 59, 11 };
  results[1] = { 57, 81, 17 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-0/test-6-polygons.off", { 0, 1 }, results);

  // Stress tests 1.

  results[0] = { 17, 17, 3 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-1/test-1-rnd-polygons-2-3.off", { 0 }, results);
  results[0] = { 16, 16, 3 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-1/test-2-rnd-polygons-2-3.off", { 0 }, results);
  results[0] = { 17, 17, 3 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-1/test-3-rnd-polygons-2-3.off", { 0 }, results);
  results[0] = { 17, 17, 3 };
  results[1] = { 19, 21, 4 };
  results[2] = { 39, 57, 13 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-1/test-4-rnd-polygons-2-4.off", { 0, 1 }, results);
  results[0] = { 12, 11, 2 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-1/test-5-rnd-polygons-1-3.off", { 0 }, results);
  results[0] = { 17, 17, 3 };
  results[1] = { 18, 20, 4 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-1/test-6-rnd-polygons-2-3.off", { 0, 1 }, results);
  results[0] = { 21, 21, 3 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-1/test-7-rnd-polygons-2-4.off", { 0 }, results);
  results[0] = { 17, 17, 3 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-1/test-8-rnd-polygons-2-10.off", { 0 }, results);
  results[0] = { 29, 33, 6 };
  results[1] = { 30, 36, 7 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-1/test-9-rnd-polygons-4-4.off", { 0, 1 }, results);
  results[0] = { 42, 53, 10 };
  results[1] = { 46, 61, 12 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-1/test-10-rnd-polygons-5-4.off", { 0, 1 }, results);

  // Stress tests 2.

  results[0] = { 17, 17, 3 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-2/test-1-rnd-polygons-2-6.off", { 0 }, results);
  results[0] = { 28, 35, 7 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-2/test-2-rnd-polygons-3-8.off", { 0 }, results);
  results[0] = { 32, 38, 7 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-2/test-3-rnd-polygons-4-4.off", { 0 }, results);
  results[0] = { 25, 27, 5 };
  results[1] = { 29, 35, 7 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-2/test-4-rnd-polygons-4-6.off", { 0, 1 }, results);
  results[0] = { 47, 61, 12 };
  results[1] = { 55, 81, 18 };
  results[2] = { 66, 100, 22 };
  results[3] = { 67, 103, 23 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-2/test-5-rnd-polygons-6-4.off", { 0, 1, 2, 3 }, results);

  results[0] = { 41, 56, 12 };
  results[1] = { 48, 69, 15 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-2/test-6-rnd-polygons-5-6.off", { 0, 1 }, results);
  results[0] = { 65, 94, 20 };
  results[1] = { 78, 123, 28 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-2/test-7-rnd-polygons-7-6.off", { 0, 1 }, results);
  results[0] = { 47, 58, 11 };
  results[1] = { 51, 66, 13 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-2/test-8-rnd-polygons-7-8.off", { 0, 1 }, results);
  results[0] = { 193, 324, 75 };
  results[1] = { 216, 389, 96 };
  results[2] = { 217, 392, 97 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-2/test-9-rnd-polygons-12-4.off", { 0, 1, 2 }, results);

  // Stress tests 3.

  results[0] = { 240, 414, 98 };
  results[1] = { 293, 540, 134 };
  results[2] = { 309, 574, 143 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-3/test-1-rnd-polygons-15-6.off", { 0, 1, 2 }, results);

  results[0] = { 407, 691, 157 };
  results[1] = { 477, 876, 214 };
  results[2] = { 489, 913, 226 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-3/test-2-rnd-polygons-20-4.off", { 0, 1, 2 }, results);

  // Stress tests 4.
  results[0] = { 401, 699, 164 };
  results[1] = { 474, 873, 214 };
  results[2] = { 488, 907, 224 };
  results[3] = { 491, 914, 226 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-4/test-1-rnd-polygons-20-6.ply", { 0, 1, 2, 3 }, results);

  results[0] = { 476, 793, 176 };
  results[1] = { 567, 1023, 245 };
  results[2] = { 602, 1100, 266 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-4/test-2-rnd-polygons-25-4.ply", { 0, 1, 2 }, results);

  results[0] = { 1275, 2521, 649 };
  results[1] = { 1445, 2988, 798 };
  results[2] = { 1501, 3113, 832 };
  results[3] = { 1519, 3155, 844 };
  run_test<Kernel, IntersectionKernel>("../data/stress-test-4/test-3-rnd-polygons-40-6.ply", { 0, 1, 2, 3 }, results);

  // Real data tests.

  results[0] = { 93, 143, 32 };
  results[1] = { 111, 185, 44 };
  results[2] = { 115, 193, 46 };
  results[3] = { 117, 197, 47 };
  run_test<Kernel, IntersectionKernel>("../data/real-data-test/test-10-polygons.off", { 0, 1, 2, 3 }, results);

  results[0] = { 222, 392, 96 };
  results[1] = { 279, 543, 143 };
  results[2] = { 290, 570, 151 };
  run_test<Kernel, IntersectionKernel>("../data/real-data-test/test-15-polygons.off", { 0, 1, 2 }, results);

  results[0] = { 370, 660, 160 };
  results[1] = { 442, 848, 217 };
  results[2] = { 446, 862, 222 };
  results[3] = { 447, 865, 223 };
  run_test<Kernel, IntersectionKernel>("../data/real-data-test/test-20-polygons.off", { 0, 1, 2, 3 }, results);

  results[0] = { 738, 1339, 327 };
  results[1] = { 919, 1760, 446 };
  results[2] = { 965, 1867, 476 };
  results[3] = { 985, 1913, 489 };
  run_test<Kernel, IntersectionKernel>("../data/real-data-test/test-40-polygons.ply", { 0, 1, 2, 3 }, results);

  timer.stop();
  std::cout << "\ntotal time: " << timer.time() << std::endl;

  const auto kernel_name = boost::typeindex::type_id<Kernel>().pretty_name();
  if (no_input != 0) {
    std::cout << kernel_name << " " << no_input << " TESTS failed due to missing input files!\n";
  }

  if (failed != 0) {
    std::cout << kernel_name << " " << failed << " TESTS failed! The space partition was empty!\n";
  }

  if (different != 0) {
    std::cout << kernel_name << " " << different << " TESTS differ from typical values!\n";
  }

  if (passed != 0)
    std::cout << kernel_name << " " << passed << " TESTS passed!\n";

  std::cout << "\n" << kernel_name << " TESTS FINISHED!" << std::endl;
}

int main(const int /* argc */, const char** /* argv */) {
  run_all_tests<SCD, EPECK>();
  //run_all_tests<EPICK, EPECK>();
  //run_all_tests<EPECK, EPECK>();
  //run_all_tests<GMPQ, GMPQ>();
  return EXIT_SUCCESS;
}
