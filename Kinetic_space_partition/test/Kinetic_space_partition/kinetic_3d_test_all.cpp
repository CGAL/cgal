#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_space_partition_3.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/IO/PLY.h>

using SCF   = CGAL::Simple_cartesian<float>;
using SCD   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
using Timer = CGAL::Real_timer;

template<typename Kernel>
bool run_test(
  const std::string input_filename,
  const std::vector<unsigned int>& ks,
  const std::vector<std::vector<unsigned int> >& results) {

  using Point_3   = typename Kernel::Point_3;
  using KSP = CGAL::Kinetic_space_partition_3<Kernel>;

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
    return false;
  }

  std::cout << input_filename << std::endl;

  for (std::size_t i = 0; i < ks.size(); i++) {
    KSP ksp(CGAL::parameters::verbose(false).debug(false));

    ksp.insert(input_vertices, input_faces);

    ksp.initialize();
    //std::cout << std::endl << "--INPUT K: " << k << std::endl;
    ksp.partition(ks[i]);

    CGAL::Linear_cell_complex_for_combinatorial_map<3, 3, CGAL::Linear_cell_complex_traits<3, CGAL::Exact_predicates_exact_constructions_kernel>, typename KSP::Linear_cell_complex_min_items> lcc;
    ksp.get_linear_cell_complex(lcc);

    std::vector<unsigned int> cells = { 0, 2, 3 }, count;
    count = lcc.count_cells(cells);

    std::cout << ksp.number_of_volumes() << std::endl;

    if (results[i][0] != count[0] || results[i][1] != count[2] || results[i][2] != count[3]) {
      std::cout << "TEST FAILED: Partitioning has not expected number of vertices, faces or volumes for k = " << ks[i] << std::endl;

      std::cout << "Expectation:" << std::endl;
      std::cout << "v: " << results[i][0] << " f : " << results[i][1] << " v : " << results[i][2] << std::endl;
      std::cout << "Result k = " << " vertices : " << count[0] << " faces : " << count[2] << " volumes : " << count[3] << std::endl;
      //assert(false);
    }
  }

  return true;
}

template<typename Kernel>
void run_all_tests() {
  std::cout.precision(10);
  std::vector< std::vector<double> > all_times;

  // All results are precomputed for k = 1!
  std::vector<std::vector<unsigned int> > results(3); //

  //run_test<Kernel>("20-inserted-polygons.ply", { 3 }, results);

  results[0] = { 58, 89, 20 };
  results[1] = { 63, 102, 24 };
  results[2] = { 63, 106, 26 };
  run_test<Kernel>("data/stress-test-5/test-2-rnd-polygons-20-4.off", { 1, 2, 3 }, results);

  results[0] = { 206, 385, 99 };
  results[1] = { 232, 449, 118 };
  results[2] = { 265, 540, 147 };
  run_test<Kernel>("data/real-data-test/test-15-polygons.off", { 1, 2, 3 }, results);

  results[0] = { 39, 49, 10 };
  results[1] = { 48, 70, 16 };
  results[2] = { 54, 84, 20 };
  run_test<Kernel>("data/edge-case-test/test-same-time.off", { 1, 2, 3 }, results);

  // Edge tests.
  results[0] = { 18, 20, 4 };
  run_test<Kernel>("data/edge-case-test/test-2-polygons.off", { 1 }, results);

  results[0] = { 22, 25, 5 };
  run_test<Kernel>("data/edge-case-test/test-4-polygons.off", { 1 }, results);

  results[0] = { 22, 25, 5 };
  run_test<Kernel>("data/edge-case-test/test-5-polygons.off", { 1 }, results);

  results[0] = { 40, 52, 11 };
  results[1] = { 51, 77, 18 };
  run_test<Kernel>("data/edge-case-test/test-local-global-1.off", { 1, 2 }, results);

  results[0] = { 40, 52, 11 };
  results[1] = { 49, 73, 17 };
  results[2] = { 54, 84, 20 };
  run_test<Kernel>("data/edge-case-test/test-local-global-2.off", { 1, 2, 3 }, results);

  // Stress tests 0.
  results[0] = { 14, 13, 2 };
  run_test<Kernel>("data/stress-test-0/test-1-polygon-a.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel>("data/stress-test-0/test-1-polygon-b.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel>("data/stress-test-0/test-1-polygon-c.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel>("data/stress-test-0/test-1-polygon-d.off", { 1 }, results);
  results[0] = { 20, 22, 4 };
  run_test<Kernel>("data/stress-test-0/test-2-polygons-ab.off", { 1 }, results);
  results[0] = { 19, 19, 3 };
  results[1] = { 20, 22, 4 };
  run_test<Kernel>("data/stress-test-0/test-2-polygons-ac.off", { 1, 2 }, results);
  results[0] = { 20, 22, 4 };
  run_test<Kernel>("data/stress-test-0/test-2-polygons-ad.off", { 1 }, results);
  results[0] = { 18, 18, 3 };
  run_test<Kernel>("data/stress-test-0/test-2-polygons-bc.off", { 1 }, results);
  results[0] = { 18, 18, 3 };
  results[1] = { 19, 21, 4 };
  run_test<Kernel>("data/stress-test-0/test-2-polygons-bd.off", { 1, 2 }, results);
  results[0] = { 19, 21, 4 };
  run_test<Kernel>("data/stress-test-0/test-2-polygons-cd.off", { 1 }, results);
  results[0] = { 27, 32, 6 };
  run_test<Kernel>("data/stress-test-0/test-3-polygons-abc.off", { 1 }, results);
  results[0] = { 28, 33, 6 };
  results[1] = { 30, 39, 8 };
  run_test<Kernel>("data/stress-test-0/test-3-polygons-abd.off", { 1, 2 }, results);
  results[0] = { 27, 32, 6 };
  results[1] = { 28, 35, 7 };
  run_test<Kernel>("data/stress-test-0/test-3-polygons-acd.off", { 1, 2 }, results);
  results[0] = { 25, 28, 5 };
  results[1] = { 26, 31, 6 };
  run_test<Kernel>("data/stress-test-0/test-3-polygons-bcd.off", { 1, 2 }, results);
  results[0] = { 36, 46, 9 };
  results[1] = { 38, 52, 11 };
  run_test<Kernel>("data/stress-test-0/test-4-polygons-abcd.off", { 1, 2 }, results);
  results[0] = { 54, 76, 16 };
  results[1] = { 64, 102, 24 };
  results[2] = { 67, 109, 26 };
  run_test<Kernel>("data/stress-test-0/test-6-polygons.off", { 1, 2, 3 }, results);

  // Stress tests 1.

  results[0] = { 14, 13, 2 };
  run_test<Kernel>("data/stress-test-1/test-1-rnd-polygons-1-4.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel>("data/stress-test-1/test-2-rnd-polygons-1-4.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel>("data/stress-test-1/test-3-rnd-polygons-1-4.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel>("data/stress-test-1/test-4-rnd-polygons-1-4.off", { 1 }, results);
  results[0] = { 18, 18, 3 };
  results[1] = { 20, 22, 4 };
  run_test<Kernel>("data/stress-test-1/test-5-rnd-polygons-2-4.off", { 1, 2 }, results);
  results[0] = { 18, 18, 3 };
  results[1] = { 19, 21, 4 };
  run_test<Kernel>("data/stress-test-1/test-6-rnd-polygons-2-4.off", { 1, 2 }, results);
  results[0] = { 18, 18, 3 };
  run_test<Kernel>("data/stress-test-1/test-7-rnd-polygons-2-4.off", { 1 }, results);

  // Stress tests 2.

  results[0] = { 14, 13, 2 };
  run_test<Kernel>("data/stress-test-2/test-1-rnd-polygons-1-4.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel>("data/stress-test-2/test-2-rnd-polygons-1-4.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel>("data/stress-test-2/test-3-rnd-polygons-1-4.off", { 1 }, results);
  results[0] = { 14, 13, 2 };
  run_test<Kernel>("data/stress-test-2/test-4-rnd-polygons-1-3.off", { 1 }, results);
  results[0] = { 17, 17, 3 };
  results[1] = { 19, 21, 4 };
  run_test<Kernel>("data/stress-test-2/test-5-rnd-polygons-2-4.off", { 1, 2 }, results);
  results[0] = { 22, 23, 4 };
  run_test<Kernel>("data/stress-test-2/test-6-rnd-polygons-3-4.off", { 1 }, results);

  // Stress tests 3.

  results[0] = { 18, 18, 3 };
  results[1] = { 20, 22, 4 };
  run_test<Kernel>("data/stress-test-3/test-1-rnd-polygons-2-3.off", { 1, 2 }, results);
  results[0] = { 17, 17, 3 };
  run_test<Kernel>("data/stress-test-3/test-2-rnd-polygons-2-3.off", { 1 }, results);
  results[0] = { 18, 18, 3 };
  results[1] = { 19, 21, 4 };
  run_test<Kernel>("data/stress-test-3/test-3-rnd-polygons-2-3.off", { 1, 2 }, results);
  results[0] = { 17, 17, 3 };
  results[1] = { 19, 21, 4 };
  run_test<Kernel>("data/stress-test-3/test-4-rnd-polygons-2-4.off", { 1, 2 }, results);
  //results[0] = { 12, 11, 2 };
  //run_test<Kernel>("data/stress-test-3/test-5-rnd-polygons-1-3.off", { 1 }, results);
  results[0] = { 18, 18, 3 };
  results[1] = { 19, 21, 4 };
  run_test<Kernel>("data/stress-test-3/test-6-rnd-polygons-2-3.off", { 1, 2 }, results);
  results[0] = { 21, 21, 3 };
  results[1] = { 22, 24, 4 };
  run_test<Kernel>("data/stress-test-3/test-7-rnd-polygons-2-4.off", { 1, 2 }, results);
  results[0] = { 17, 17, 3 };
  results[1] = { 18, 20, 4 };
  run_test<Kernel>("data/stress-test-3/test-8-rnd-polygons-2-10.off", { 1, 2 }, results);
  results[0] = { 31, 37, 7 };
  results[1] = { 34, 46, 10 };
  results[2] = { 39, 57, 13 };
  run_test<Kernel>("data/stress-test-3/test-9-rnd-polygons-4-4.off", { 1, 2, 3 }, results);
  results[0] = { 55, 84, 19 };
  run_test<Kernel>("data/stress-test-3/test-10-rnd-polygons-5-4.off", { 1 }, results);

  // Stress tests 4.

  results[0] = { 17, 17, 3 };
  run_test<Kernel>("data/stress-test-4/test-1-rnd-polygons-2-6.off", { 1 }, results);
  results[0] = { 27, 32, 6 };
  results[1] = { 29, 38, 8 };
  run_test<Kernel>("data/stress-test-4/test-2-rnd-polygons-3-8.off", { 1, 2 }, results);
  results[0] = { 32, 38, 7 };
  results[1] = { 35, 45, 9 };
  results[2] = { 37, 51, 11 };
  run_test<Kernel>("data/stress-test-4/test-3-rnd-polygons-4-4.off", { 1, 2, 3 }, results);
  results[0] = { 25, 27, 5 };
  results[1] = { 33, 41, 8 };
  results[2] = { 37, 53, 12 };
  run_test<Kernel>("data/stress-test-4/test-4-rnd-polygons-4-6.off", { 1, 2, 3 }, results);
  results[0] = { 61, 91, 20 };
  results[1] = { 73, 121, 29 };
  results[2] = { 83, 145, 36 };
  run_test<Kernel>("data/stress-test-4/test-5-rnd-polygons-6-4.off", { 1, 2, 3 }, results);

  results[0] = { 45, 62, 13 };
  results[1] = { 50, 75, 17 };
  run_test<Kernel>("data/stress-test-4/test-6-rnd-polygons-5-6.off", { 1, 2 }, results);
  results[0] = { 64, 97, 22 };
  results[1] = { 84, 141, 34 };
  results[2] = { 88, 151, 37 };
  run_test<Kernel>("data/stress-test-4/test-7-rnd-polygons-7-6.off", { 1, 2, 3 }, results);
  results[0] = { 56, 77, 16 };
  results[1] = { 68, 107, 25 };
  results[2] = { 69, 110, 26 };
  run_test<Kernel>("data/stress-test-4/test-8-rnd-polygons-7-8.off", { 1, 2, 3 }, results);
  results[0] = { 172, 304, 74 };
  results[1] = { 192, 366, 95 };
  results[2] = { 198, 382, 100};
  run_test<Kernel>("data/stress-test-4/test-9-rnd-polygons-12-4.off", { 1, 2, 3 }, results);

  // Stress tests 5.

  results[0] = { 202, 351, 84 };
  results[1] = { 232, 427, 107 };
  results[2] = { 284, 558, 146 };
  run_test<Kernel>("data/stress-test-5/test-1-rnd-polygons-15-6.off", { 1, 2, 3 }, results);
  results[0] = { 58, 89, 20 };
  results[1] = { 63, 102, 24 };
  results[2] = { 63, 106, 26 };
  run_test<Kernel>("data/stress-test-5/test-2-rnd-polygons-20-4.off", { 1, 2, 3 }, results);

  // Real data tests.

  results[0] = { 91, 143, 33 };
  results[1] = { 109, 187, 46 };
  results[2] = { 127, 233, 60 };
  run_test<Kernel>("data/real-data-test/test-10-polygons.off", { 1, 2, 3 }, results);

  results[0] = { 1006, 2067, 552 };
  results[1] = { 973, 1984, 527 };
  results[2] = { 1186, 2560, 708 };
  run_test<Kernel>("data/real-data-test/test-40-polygons.ply", { 1, 2, 3 }, results);

  const auto kernel_name = boost::typeindex::type_id<Kernel>().pretty_name();
  std::cout << std::endl << kernel_name << " TESTS SUCCESS!" << std::endl << std::endl;
}

#include <CGAL/intersections.h>

int main(const int /* argc */, const char** /* argv */) {
  // run_all_tests<SCF>();
  // run_all_tests<SCD>();
  //run_all_tests<EPECK>();

  // Passes all tests except for those when
  // intersections lead to accumulated errors.
  //build();
  run_all_tests<EPICK>();
  return EXIT_SUCCESS;
}
