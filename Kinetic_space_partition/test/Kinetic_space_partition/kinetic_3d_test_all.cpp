#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/Kinetic_space_partition_3.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/IO/PLY.h>

using SCD   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

std::size_t different = 0;

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
    different++;
    return false;
  }

  std::cout << input_filename << std::endl;

  for (std::size_t i = 0; i < ks.size(); i++) {
    KSP ksp(CGAL::parameters::verbose(false).debug(false));

    ksp.insert(input_vertices, input_faces);

    ksp.initialize();
    //std::cout << std::endl << "--INPUT K: " << k << std::endl;
    ksp.partition(ks[i]);

    CGAL::Linear_cell_complex_for_combinatorial_map<3, 3, CGAL::Linear_cell_complex_traits<3, Kernel>, typename KSP::Linear_cell_complex_min_items> lcc;
    ksp.get_linear_cell_complex(lcc);

    std::vector<unsigned int> cells = { 0, 2, 3 }, count;
    count = lcc.count_cells(cells);

    std::cout << ksp.number_of_volumes() << std::endl;

    if (results[i][0] != count[0] || results[i][1] != count[2] || results[i][2] != count[3]) {
      std::cout << "TEST differs: Partitioning has not expected number of vertices, faces or volumes for k = " << ks[i] << std::endl;

      std::cout << "Expectation:" << std::endl;
      std::cout << "v: " << results[i][0] << " f : " << results[i][1] << " v : " << results[i][2] << std::endl;
      std::cout << "Result k = " << " vertices : " << count[0] << " faces : " << count[2] << " volumes : " << count[3] << std::endl;
      std::cout << input_filename << std::endl;

      different++;
    }
    else std::cout << "TEST PASSED k = " << ks[i] << " " << input_filename << std::endl;
  }

  return true;
}

template<typename Kernel, typename IntersectionKernel>
void run_all_tests() {
  different = 0;
  std::cout.precision(10);
  std::vector< std::vector<double> > all_times;

  std::vector<std::vector<unsigned int> > results(3);

  results[0] = { 40, 52, 11 };
  results[1] = { 48, 70, 16 };
  results[2] = { 54, 84, 20 };
  run_test<Kernel, IntersectionKernel>("data/edge-case-test/test-same-time.off", { 1, 2, 3 }, results);

  // Edge tests.
  results[0] = { 18, 20, 4 };
  run_test<Kernel, IntersectionKernel>("data/edge-case-test/test-2-polygons.off", { 1 }, results);

  results[0] = { 22, 25, 5 };
  run_test<Kernel, IntersectionKernel>("data/edge-case-test/test-4-polygons.off", { 1 }, results);

  results[0] = { 22, 25, 5 };
  run_test<Kernel, IntersectionKernel>("data/edge-case-test/test-5-polygons.off", { 1 }, results);

  results[0] = { 39, 49, 10 };
  results[1] = { 49, 73, 17 };
  run_test<Kernel, IntersectionKernel>("data/edge-case-test/test-local-global-1.off", { 1, 2 }, results);

  results[0] = { 39, 49, 10 };
  results[1] = { 51, 77, 18 };
  results[2] = { 54, 84, 20 };
  run_test<Kernel, IntersectionKernel>("data/edge-case-test/test-local-global-2.off", { 1, 2, 3 }, results);

  // Stress tests 0.
  results[0] = { 14, 13, 2 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-1-polygon-a.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-1-polygon-b.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-1-polygon-c.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-1-polygon-d.off", { 1 }, results);
  results[0] = { 18, 18, 3 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-2-polygons-ab.off", { 1 }, results);
  results[0] = { 19, 19, 3 };
  results[1] = { 20, 22, 4 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-2-polygons-ac.off", { 1, 2 }, results);
  results[0] = { 19, 19, 3 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-2-polygons-ad.off", { 1 }, results);
  results[0] = { 18, 18, 3 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-2-polygons-bc.off", { 1 }, results);
  results[0] = { 18, 18, 3 };
  results[1] = { 19, 21, 4 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-2-polygons-bd.off", { 1, 2 }, results);
  results[0] = { 19, 21, 4 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-2-polygons-cd.off", { 1 }, results);
  results[0] = { 26, 29, 5 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-3-polygons-abc.off", { 1 }, results);
  results[0] = { 28, 31, 5 };
  results[1] = { 30, 39, 8 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-3-polygons-abd.off", { 1, 2 }, results);
  results[0] = { 25, 28, 5 };
  results[1] = { 27, 32, 6 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-3-polygons-acd.off", { 1, 2 }, results);
  results[0] = { 25, 28, 5 };
  results[1] = { 26, 31, 6 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-3-polygons-bcd.off", { 1, 2 }, results);
  results[0] = { 34, 42, 8 };
  results[1] = { 38, 52, 11 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-4-polygons-abcd.off", { 1, 2 }, results);
  results[0] = { 50, 68, 14 };
  results[1] = { 56, 82, 18 };
  results[2] = { 67, 109, 26 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-0/test-6-polygons.off", { 1, 2, 3 }, results);

 // Stress tests 1.

  results[0] = { 14, 13, 2 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-1/test-1-rnd-polygons-1-4.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-1/test-2-rnd-polygons-1-4.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-1/test-3-rnd-polygons-1-4.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-1/test-4-rnd-polygons-1-4.off", { 1 }, results);
  results[0] = { 18, 18, 3 };
  results[1] = { 20, 22, 4 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-1/test-5-rnd-polygons-2-4.off", { 1, 2 }, results);
  results[0] = { 18, 18, 3 };
  results[1] = { 19, 21, 4 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-1/test-6-rnd-polygons-2-4.off", { 1, 2 }, results);
  results[0] = { 18, 18, 3 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-1/test-7-rnd-polygons-2-4.off", { 1 }, results);

  // Stress tests 2.

  results[0] = { 14, 13, 2 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-2/test-1-rnd-polygons-1-4.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-2/test-2-rnd-polygons-1-4.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-2/test-3-rnd-polygons-1-4.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-2/test-4-rnd-polygons-1-3.off", { 1 }, results);
  results[0] = { 17, 17, 3 };
  results[1] = { 19, 21, 4 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-2/test-5-rnd-polygons-2-4.off", { 1, 2 }, results);
  results[0] = { 24, 27, 5 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-2/test-6-rnd-polygons-3-4.off", { 1 }, results);

  // Stress tests 3.

  results[0] = { 18, 18, 3 };
  results[1] = { 20, 22, 4 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-3/test-1-rnd-polygons-2-3.off", { 1, 2 }, results);
  results[0] = { 17, 17, 3 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-3/test-2-rnd-polygons-2-3.off", { 1 }, results);
  results[0] = { 18, 18, 3 };
  results[1] = { 19, 21, 4 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-3/test-3-rnd-polygons-2-3.off", { 1, 2 }, results);
  results[0] = { 17, 17, 3 };
  results[1] = { 19, 21, 4 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-3/test-4-rnd-polygons-2-4.off", { 1, 2 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-3/test-5-rnd-polygons-1-3.off", { 1 }, results);
  results[0] = { 18, 18, 3 };
  results[1] = { 19, 21, 4 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-3/test-6-rnd-polygons-2-3.off", { 1, 2 }, results);
  results[0] = { 21, 21, 3 };
  results[1] = { 22, 24, 4 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-3/test-7-rnd-polygons-2-4.off", { 1, 2 }, results);
  results[0] = { 17, 17, 3 };
  results[1] = { 18, 20, 4 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-3/test-8-rnd-polygons-2-10.off", { 1, 2 }, results);
  results[0] = { 31, 37, 7 };
  results[1] = { 34, 46, 10 };
  results[2] = { 39, 57, 13 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-3/test-9-rnd-polygons-4-4.off", { 1, 2, 3 }, results);
  results[0] = { 51, 72, 15 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-3/test-10-rnd-polygons-5-4.off", { 1 }, results);

  // Stress tests 4.

  results[0] = { 17, 17, 3 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-4/test-1-rnd-polygons-2-6.off", { 1 }, results);
  results[0] = { 27, 32, 6 };
  results[1] = { 29, 38, 8 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-4/test-2-rnd-polygons-3-8.off", { 1, 2 }, results);
  results[0] = { 32, 38, 7 };
  results[1] = { 35, 45, 9 };
  results[2] = { 37, 51, 11 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-4/test-3-rnd-polygons-4-4.off", { 1, 2, 3 }, results);
  results[0] = { 25, 27, 5 };
  results[1] = { 30, 36, 7 };
  results[2] = { 37, 53, 12 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-4/test-4-rnd-polygons-4-6.off", { 1, 2, 3 }, results);
  results[0] = { 61, 91, 20 };
  results[1] = { 73, 121, 29 };
  results[2] = { 83, 145, 36 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-4/test-5-rnd-polygons-6-4.off", { 1, 2, 3 }, results);

  results[0] = { 43, 58, 12 };
  results[1] = { 50, 75, 17 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-4/test-6-rnd-polygons-5-6.off", { 1, 2 }, results);
  results[0] = { 63, 94, 21 };
  results[1] = { 84, 141, 34 };
  results[2] = { 90, 157, 39 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-4/test-7-rnd-polygons-7-6.off", { 1, 2, 3 }, results);
  results[0] = { 56, 77, 16 };
  results[1] = { 68, 107, 25 };
  results[2] = { 69, 110, 26 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-4/test-8-rnd-polygons-7-8.off", { 1, 2, 3 }, results);
  results[0] = { 173, 305, 74 };
  results[1] = { 194, 370, 96 };
  results[2] = { 207, 407, 108 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-4/test-9-rnd-polygons-12-4.off", { 1, 2, 3 }, results);

  // Stress tests 5.

  results[0] = { 185, 308, 71 };
  results[1] = { 223, 406, 101 };
  results[2] = { 277, 548, 145 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-5/test-1-rnd-polygons-15-6.off", { 1, 2, 3 }, results);

  results[0] = { 50, 71, 15 };
  results[1] = { 56, 85, 19 };
  results[2] = { 63, 102, 24 };
  run_test<Kernel, IntersectionKernel>("data/stress-test-5/test-2-rnd-polygons-20-4.off", { 1, 2, 3 }, results);

  // Real data tests.

  results[0] = { 94, 150, 35 };
  results[1] = { 111, 191, 47 };
  results[2] = { 127, 233, 60 };
  run_test<Kernel, IntersectionKernel>("data/real-data-test/test-10-polygons.off", { 1, 2, 3 }, results);

  results[0] = { 206, 385, 99 };
  results[1] = { 237, 462, 122 };
  results[2] = { 260, 529, 144 };
  run_test<Kernel, IntersectionKernel>("data/real-data-test/test-15-polygons.off", { 1, 2, 3 }, results);

  results[0] = { 1156, 2466, 677 };
  results[1] = { 1131, 2387, 650 };
  results[2] = { 1395, 3115, 882 };
  run_test<Kernel, IntersectionKernel>("data/real-data-test/test-40-polygons.ply", { 1, 2, 3 }, results);

  const auto kernel_name = boost::typeindex::type_id<Kernel>().pretty_name();
  if (different != 0) {
    std::cout << std::endl << kernel_name << " " << different << " TESTS differ from typical values!" << std::endl << std::endl;
  }
  std::cout << std::endl << kernel_name << " TESTS SUCCESS!" << std::endl << std::endl;
}

int main(const int /* argc */, const char** /* argv */) {
  run_all_tests<SCD, EPECK>();
  run_all_tests<EPICK, EPECK>();
  //run_all_tests<EPECK, EPECK>();
  //run_all_tests<GMPQ, GMPQ>();
  return EXIT_SUCCESS;
}
