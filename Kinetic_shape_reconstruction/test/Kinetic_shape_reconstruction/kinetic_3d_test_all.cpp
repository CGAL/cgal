#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_reconstruction_3.h>
#include <CGAL/IO/OFF_reader.h>

using EPICK   = CGAL::Exact_predicates_inexact_constructions_kernel;
using Kernel  = EPICK;
using Point_3 = typename Kernel::Point_3;
using KSR     = CGAL::Kinetic_shape_reconstruction_3<Kernel>;

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

const bool run_test(
  const std::string input_filename,
  const std::vector<unsigned int>& ks,
  const std::size_t num_iters,
  std::size_t& num_tests) {

  ++num_tests;
  std::ifstream input_file(input_filename);
  std::vector<Point_3> input_vertices;
  std::vector< std::vector<std::size_t> > input_faces;
  assert(CGAL::read_OFF(input_file, input_vertices, input_faces));

  std::cout << std::endl;
  std::cout << "--INPUT FILE: " << input_filename << std::endl;
  const Polygon_map polygon_map(input_vertices);
  for (const auto k : ks) {
    std::cout << std::endl << "--INPUT K: " << k << std::endl;
    for (std::size_t iter = 0; iter < num_iters; ++iter) {
      std::cout << std::endl << "--ITERATION #" << iter + 1 << " BEGIN!" << std::endl;
      KSR ksr(false, false);
      assert(ksr.partition(input_faces, polygon_map, k));
      ksr.clear();
      std::cout << std::endl << "--ITERATION #" << iter + 1 << " END!" << std::endl;
    }
  }
  return true;
}

int main (const int argc, const char** argv) {

  std::size_t num_tests = 0;
  const std::size_t num_iters = 3;

  std::vector<unsigned int> ks;
  for (unsigned int k = 1; k <= 6; ++k) {
    ks.push_back(k);
  }
  ks.push_back(100);

  // Stress tests 0.
  assert(run_test("data/stress-test-0/test-1-polygon-a.off"    , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-1-polygon-b.off"    , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-1-polygon-c.off"    , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-1-polygon-d.off"    , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-2-polygons-ab.off"  , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-2-polygons-ac.off"  , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-2-polygons-ad.off"  , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-2-polygons-bc.off"  , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-2-polygons-bd.off"  , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-2-polygons-cd.off"  , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-3-polygons-abc.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-3-polygons-abd.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-3-polygons-acd.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-3-polygons-bcd.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-4-polygons-abcd.off", ks, num_iters, num_tests));
  assert(run_test("data/stress-test-0/test-6-polygons.off"     , ks, num_iters, num_tests));

  // Stress tests 1.
  assert(run_test("data/stress-test-1/test-1-rnd-polygons-1-4.off", ks, num_iters, num_tests));
  assert(run_test("data/stress-test-1/test-2-rnd-polygons-1-4.off", ks, num_iters, num_tests));
  assert(run_test("data/stress-test-1/test-3-rnd-polygons-1-4.off", ks, num_iters, num_tests));
  assert(run_test("data/stress-test-1/test-4-rnd-polygons-1-4.off", ks, num_iters, num_tests));
  assert(run_test("data/stress-test-1/test-5-rnd-polygons-2-4.off", ks, num_iters, num_tests));
  assert(run_test("data/stress-test-1/test-6-rnd-polygons-2-4.off", ks, num_iters, num_tests));
  assert(run_test("data/stress-test-1/test-7-rnd-polygons-2-4.off", ks, num_iters, num_tests));
  assert(run_test("data/stress-test-1/test-8-rnd-polygons-3-4.off", ks, num_iters, num_tests));

  // Stress tests 2.
  assert(run_test("data/stress-test-2/test-1-rnd-polygons-1-4.off", ks, num_iters, num_tests));
  assert(run_test("data/stress-test-2/test-2-rnd-polygons-1-4.off", ks, num_iters, num_tests));
  assert(run_test("data/stress-test-2/test-3-rnd-polygons-1-4.off", ks, num_iters, num_tests));
  assert(run_test("data/stress-test-2/test-4-rnd-polygons-1-3.off", ks, num_iters, num_tests));
  assert(run_test("data/stress-test-2/test-5-rnd-polygons-2-4.off", ks, num_iters, num_tests));
  assert(run_test("data/stress-test-2/test-6-rnd-polygons-3-4.off", ks, num_iters, num_tests));

  // Stress tests 3.
  assert(run_test("data/stress-test-3/test-1-rnd-polygons-2-3.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-2-rnd-polygons-2-3.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-3-rnd-polygons-2-3.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-4-rnd-polygons-2-4.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-5-rnd-polygons-1-3.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-6-rnd-polygons-2-3.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-7-rnd-polygons-2-4.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-8-rnd-polygons-2-10.off", ks, num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-9-rnd-polygons-4-4.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-3/test-10-rnd-polygons-5-4.off", ks, num_iters, num_tests));

  // Stress tests 4.
  assert(run_test("data/stress-test-4/test-1-rnd-polygons-2-6.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-2-rnd-polygons-3-8.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-3-rnd-polygons-4-4.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-4-rnd-polygons-4-6.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-5-rnd-polygons-6-4.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-6-rnd-polygons-5-6.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-7-rnd-polygons-7-6.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-8-rnd-polygons-7-8.off" , ks, num_iters, num_tests));
  assert(run_test("data/stress-test-4/test-9-rnd-polygons-12-4.off", ks, num_iters, num_tests));

  // Real data tests.
  assert(run_test("data/real-data-test/building-a-10polygons-10planes.off", ks, num_iters, num_tests));
  assert(run_test("data/real-data-test/building-b-15squares-15planes.off" , ks, num_iters, num_tests));

  // Edge case tests.
  assert(run_test("data/edge-case-test/test-2-polygons.off" , ks, num_iters, num_tests)); // edge touch
  assert(run_test("data/edge-case-test/test-4-polygons.off" , ks, num_iters, num_tests)); // edge touch / 2 coplanar
  assert(run_test("data/edge-case-test/test-5-polygons.off" , ks, num_iters, num_tests)); // edge touch / vertex touch / 2 coplanar

  // std::vector<unsigned int> ts;
  // ts.push_back(1); ts.push_back(2); ts.push_back(4);
  // ts.push_back(5); ts.push_back(6); ts.push_back(100);
  // assert(run_test("data/edge-case-test/test-20-polygons.off", ts, num_iters, num_tests)); // 2 overlap and coplanar

  // std::cout << std::endl << "--OUTPUT STATS:" << std::endl;
  // std::cout << "* number of iterations per test: " << num_iters << std::endl;
  // std::cout << "* k intersections: {";
  // for (const auto k : ks) {
  //   std::cout << k << ",";
  // }
  // std::cout << "...}" << std::endl;

  std::cout << std::endl << "ALL " << num_tests << " TESTS SUCCESS!" << std::endl;
  return EXIT_SUCCESS;
}
