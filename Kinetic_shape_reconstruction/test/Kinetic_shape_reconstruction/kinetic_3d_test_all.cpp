#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic_shape_partitioning_3.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/IO/PLY.h>

using SCF   = CGAL::Simple_cartesian<float>;
using SCD   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
using Timer = CGAL::Real_timer;

using Traits = typename CGAL::Kinetic_shape_partitioning_Traits_3<EPICK, EPECK, std::vector<typename EPICK::Point_3>, CGAL::Identity_property_map<typename EPICK::Point_3> >;

template<typename Point>
struct Polygon_map {

  using key_type   = std::vector<std::size_t>;
  using value_type = std::vector<Point>;
  using reference  = value_type;
  using category   = boost::readable_property_map_tag;

  const std::vector<Point>& points;
  Polygon_map(
    const std::vector<Point>& vertices) :
  points(vertices)
  { }

  friend reference get(const Polygon_map& map, const key_type& face) {
    reference polygon;
    polygon.reserve(face.size());
    std::transform(
      face.begin(), face.end(),
      std::back_inserter(polygon),
      [&](const std::size_t vertex_index) -> Point {
        return map.points[vertex_index];
      });
    return polygon;
  }
};

template<typename Traits>
bool run_test(
  const std::string input_filename,
  const std::vector<unsigned int>& ks,
  const std::size_t num_iters,
  const std::vector<int>& results,
  std::vector< std::vector<double> >& all_times,
  std::size_t& num_tests) {

  using Point_3   = typename Traits::Kernel::Point_3;
  using Segment_3 = typename Traits::Kernel::Segment_3;

  using Surface_mesh = CGAL::Surface_mesh<Point_3>;
  using KSP = CGAL::Kinetic_shape_partitioning_3<Traits>;

  ++num_tests;
  std::string baseDir = "C:/dev/kinetic/Kinetic_shape_reconstruction/examples/Kinetic_shape_reconstruction/";
  std::string filename = baseDir + input_filename;
  std::ifstream input_file_off(filename);
  std::ifstream input_file_ply(filename);
  std::vector<Point_3> input_vertices;
  std::vector< std::vector<std::size_t> > input_faces;

  if (CGAL::IO::read_OFF(input_file_off, input_vertices, input_faces)) {
    std::cout << "* reading the OFF file: " << filename << "!" << std::endl;
    input_file_off.close();
  } else if (CGAL::IO::read_PLY(input_file_ply, input_vertices, input_faces)) {
    std::cout << "* reading the PLY file: " << filename << "!" << std::endl;
    input_file_ply.close();
  } else {
    std::cerr << "ERROR: can't read the OFF/PLY file " << filename << "!" << std::endl;
    return false;
  }
  std::vector<double> times;

  std::cout << std::endl;
  std::cout << "--INPUT FILE: " << filename << std::endl;
  const Polygon_map<Point_3> polygon_map(input_vertices);
  for (const unsigned int k : ks) {
    std::cout << std::endl << "--INPUT K: " << k << std::endl;

    double time = 0.0;
    for (std::size_t iter = 0; iter < num_iters; ++iter) {
      std::cout << std::endl << "--ITERATION #" << iter + 1 << " BEGIN!" << std::endl;
       KSP ksp(true, false); // first verbose, second debug

      // Running KSR.
      Timer timer;
      timer.start();

      bool is_ksp_success = ksp.initialize(
        input_faces, polygon_map);

      if (is_ksp_success)
        ksp.partition(k);

      assert(is_ksp_success);
      if (!is_ksp_success) return false;
      timer.stop();
      time += timer.time();

      // Testing results.

      const int num_support_planes = ksp.number_of_support_planes();

      const int num_vertices = static_cast<int>(ksp.number_of_vertices());
      const int num_faces = static_cast<int>(ksp.number_of_faces());
      const int num_volumes  = static_cast<int>(ksp.number_of_volumes());

      std::cout << std::endl << "--RESULTS: ";
      std::cout << num_support_planes << ",";

      std::cout << num_vertices << ",";
      std::cout << num_faces    << ",";
      std::cout << num_volumes  << std::endl;

/*
      assert(num_support_planes > 6);

      if (num_support_planes <= 6) return false;

      assert(results.size() == 6);
      assert(num_support_planes == results[0]);

      if (results.size() != 6) return false;
      if (num_support_planes != results[0]) return false;

      assert(num_vertices == results[2]);
      assert(num_edges    == results[3]);
      assert(num_faces    >= results[4]);
      assert(num_volumes  >= results[5]);

      if (num_vertices != results[2]) return false;
      if (num_edges    != results[3]) return false;
      if (num_faces     < results[4]) return false;
      if (num_volumes   < results[5]) return false;*/

      CGAL::Linear_cell_complex_for_combinatorial_map<3, 3> lcc;
      ksp.get_linear_cell_complex(lcc);
/*
      std::vector<Point_3> output_vertices;
      ksp.output_partition_vertices(
        std::back_inserter(output_vertices));
      assert(static_cast<std::size_t>(num_vertices) == output_vertices.size());
      if (static_cast<std::size_t>(num_vertices) != output_vertices.size()) return false;


      std::vector<Segment_3> output_edges;
      ksr.output_partition_edges(
        std::back_inserter(output_edges));
      assert(static_cast<std::size_t>(num_edges) == output_edges.size());
      if (static_cast<std::size_t>(num_edges) != output_edges.size()) return false;

      std::vector< std::vector<std::size_t> > output_faces;
      ksr.output_partition_faces(
        std::back_inserter(output_faces));
      assert(static_cast<std::size_t>(num_faces) == output_faces.size());
      if (static_cast<std::size_t>(num_faces) != output_faces.size()) return false;

      std::vector<Surface_mesh> output_volumes;
      ksr.output_partition_volumes(
        std::back_inserter(output_volumes));
      assert(static_cast<std::size_t>(num_volumes) == output_volumes.size());
      if (static_cast<std::size_t>(num_volumes) != output_volumes.size()) return false;*/

      ksp.clear();
      assert(ksp.number_of_support_planes() == 0);
      assert(ksp.number_of_vertices()       == 0);
      assert(ksp.number_of_faces()          == 0);
      assert(ksp.number_of_volumes()        == 0);

      if (ksp.number_of_support_planes() != 0) return false;
      if (ksp.number_of_vertices()       != 0) return false;
      if (ksp.number_of_faces()          != 0) return false;
      if (ksp.number_of_volumes()        != 0) return false;

      std::cout << std::endl << "--ITERATION #" << iter + 1 << " END!" << std::endl;
    }
    time /= static_cast<double>(num_iters);
    times.push_back(time);
  }

  assert(times.size() == ks.size());
  if (times.size() != ks.size()) return false;
  all_times.push_back(times);
  return true;
}

template<typename Traits>
void run_all_tests() {
  std::size_t num_tests = 0;
  const std::size_t num_iters = 1;

  std::cout.precision(10);
  std::vector< std::vector<double> > all_times;

  // All results are precomputed for k = 1!
  std::vector<int> results;

  // Number of allowed intersections k.
  std::vector<unsigned int> ks;
  for (unsigned int k = 1; k <= 6; ++k) {
    ks.push_back(k);
  }
  results = { 9,1,28,56,35,6 };
  run_test<Traits>("data/stress-test-1/test-8-rnd-polygons-3-4.off", ks, num_iters, results, all_times, num_tests);
  results = { 16,1,133,315,212,34 };
  run_test<Traits>("data/real-data-test/test-10-polygons.off", ks, num_iters, results, all_times, num_tests);
  results = { 10,1,37,77,46,6 };
  run_test<Traits>("data/stress-test-4/test-4-rnd-polygons-4-6.off", ks, num_iters, results, all_times, num_tests);
  results = { 10,1,37,77,46,6 };
  run_test<Traits>("data/edge-case-test/test-box.off", ks, num_iters, results, all_times, num_tests);
  results = {7,1,12,20,11,2};
  run_test<Traits>("data/edge-case-test/test-flat-bbox-xy-split.off", ks, num_iters, results, all_times, num_tests);

  // Edge case tests.

  // flat bbox / 2 coplanar in XY
  //results = { 7,1,14,24,13,2 };
  //assert(run_test<Traits>("data/stress-test-0/test-1-polygon-a.off", ks, num_iters, results, all_times, num_tests));

  //results = { 8,1,20,37,21,3 };
  //assert(run_test<Traits>("data/stress-test-0/test-2-polygons-ab.off", ks, num_iters, results, all_times, num_tests));
  //results = {7,1,12,20,11,2};
  //assert(run_test<Traits>("data/edge-case-test/test-flat-bbox-xy-split.off", ks, num_iters, results, all_times, num_tests));
  //results = { 10,1,38,78,46,6 };
  //assert(run_test<Traits>("data/stress-test-0/test-4-polygons-abcd.off", ks, num_iters, results, all_times, num_tests));


  // flat bbox / 2 coplanar in XZ
 /* results = {7,1,12,20,11,2};
  assert(run_test<Traits>("data/edge-case-test/test-flat-bbox-xz.off", ks, num_iters, results, all_times, num_tests));

  // flat bbox / 2 coplanar in YZ
  results = {7,1,12,20,11,2};
  assert(run_test<Traits>("data/edge-case-test/test-flat-bbox-yz.off", ks, num_iters, results, all_times, num_tests));

  // edge touch
  results = {8,1,18,33,19,3};
  assert(run_test<Traits>("data/edge-case-test/test-2-polygons.off"  , ks, num_iters, results, all_times, num_tests));

  // edge touch / 2 coplanar
  results = {9,1,24,46,27,4};
  assert(run_test<Traits>("data/edge-case-test/test-4-polygons.off"  , ks, num_iters, results, all_times, num_tests));

  // edge touch / vertex touch / 2 coplanar
  results = {9,1,24,46,27,4};
  assert(run_test<Traits>("data/edge-case-test/test-5-polygons.off"  , ks, num_iters, results, all_times, num_tests));

  // multiple collinear and duplicate input vertices
  results = {8,1,18,33,19,3};
  assert(run_test<Traits>("data/edge-case-test/test-collinear.off"     , ks, num_iters, results, all_times, num_tests));

  // all events happen at the same time
  results = {12,1,54,117,74,11};
  assert(run_test<Traits>("data/edge-case-test/test-same-time.off"     , ks, num_iters, results, all_times, num_tests));

  // failure case #1 that produces holes
  results = {12,1,54,117,69,9};
  assert(run_test<Traits>("data/edge-case-test/test-local-global-1.off", ks, num_iters, results, all_times, num_tests));

  // failure case #2 that produces holes
  results = {12,1,54,117,70,9};
  assert(run_test<Traits>("data/edge-case-test/test-local-global-2.off", ks, num_iters, results, all_times, num_tests));*/

  // Stress tests 0.
  //results = {7,1,14,24,13,2};
  //assert(run_test<Traits>("data/stress-test-0/test-1-polygon-a.off"    , ks, num_iters, results, all_times, num_tests));
  /*results = {7,1,14,24,13,2};
  assert(run_test<Traits>("data/stress-test-0/test-1-polygon-b.off"    , ks, num_iters, results, all_times, num_tests));
  results = {7,1,14,24,13,2};
  assert(run_test<Traits>("data/stress-test-0/test-1-polygon-c.off"    , ks, num_iters, results, all_times, num_tests));
  results = {7,1,14,24,13,2};
  assert(run_test<Traits>("data/stress-test-0/test-1-polygon-d.off"    , ks, num_iters, results, all_times, num_tests));
  results = {8,1,20,37,21,3};
  assert(run_test<Traits>("data/stress-test-0/test-2-polygons-ab.off"  , ks, num_iters, results, all_times, num_tests));
  results = {8,1,20,37,21,3};
  assert(run_test<Traits>("data/stress-test-0/test-2-polygons-ac.off"  , ks, num_iters, results, all_times, num_tests));
  results = {8,1,20,37,21,3};
  assert(run_test<Traits>("data/stress-test-0/test-2-polygons-ad.off"  , ks, num_iters, results, all_times, num_tests));
  results = {8,1,18,32,18,3};
  assert(run_test<Traits>("data/stress-test-0/test-2-polygons-bc.off"  , ks, num_iters, results, all_times, num_tests));
  results = {8,1,19,35,20,3};
  assert(run_test<Traits>("data/stress-test-0/test-2-polygons-bd.off"  , ks, num_iters, results, all_times, num_tests));
  results = {8,1,19,35,21,4};
  assert(run_test<Traits>("data/stress-test-0/test-2-polygons-cd.off"  , ks, num_iters, results, all_times, num_tests));
  results = {9,1,27,52,30,4};
  assert(run_test<Traits>("data/stress-test-0/test-3-polygons-abc.off" , ks, num_iters, results, all_times, num_tests));

  results = {9,1,30,60,34,4};
  assert(run_test<Traits>("data/stress-test-0/test-3-polygons-abd.off" , ks, num_iters, results, all_times, num_tests));
  results = {9,1,28,55,33,5};
  assert(run_test<Traits>("data/stress-test-0/test-3-polygons-acd.off" , ks, num_iters, results, all_times, num_tests));*/
  /*results = {9,1,26,50,30,5};
  assert(run_test<Traits>("data/stress-test-0/test-3-polygons-bcd.off" , ks, num_iters, results, all_times, num_tests));
  results = {10,1,38,78,46,6};
  assert(run_test<Traits>("data/stress-test-0/test-4-polygons-abcd.off", ks, num_iters, results, all_times, num_tests));
  results = { 12,1,67,149,90,11 };
  assert(run_test<Traits>("data/stress-test-0/test-6-polygons.off", ks, num_iters, results, all_times, num_tests));*/

  // Stress tests 1.
/*
  results = {7,1,14,24,13,2};
  assert(run_test<Traits>("data/stress-test-1/test-1-rnd-polygons-1-4.off", ks, num_iters, results, all_times, num_tests));
  results = {7,1,14,24,13,2};
  assert(run_test<Traits>("data/stress-test-1/test-2-rnd-polygons-1-4.off", ks, num_iters, results, all_times, num_tests));
  results = {7,1,14,24,13,2};
  assert(run_test<Traits>("data/stress-test-1/test-3-rnd-polygons-1-4.off", ks, num_iters, results, all_times, num_tests));
  results = {7,1,14,24,13,2};
  assert(run_test<Traits>("data/stress-test-1/test-4-rnd-polygons-1-4.off", ks, num_iters, results, all_times, num_tests));
  results = {8,1,20,37,21,3};
  assert(run_test<Traits>("data/stress-test-1/test-5-rnd-polygons-2-4.off", ks, num_iters, results, all_times, num_tests));
  results = {8,1,19,35,20,3};
  assert(run_test<Traits>("data/stress-test-1/test-6-rnd-polygons-2-4.off", ks, num_iters, results, all_times, num_tests));
  results = {8,1,20,37,22,4};
  assert(run_test<Traits>("data/stress-test-1/test-7-rnd-polygons-2-4.off", ks, num_iters, results, all_times, num_tests));*/

  // Stress tests 2.
/*
  results = {7,1,14,24,13,2};
  assert(run_test<Traits>("data/stress-test-2/test-1-rnd-polygons-1-4.off", ks, num_iters, results, all_times, num_tests));
  results = {7,1,14,24,13,2};
  assert(run_test<Traits>("data/stress-test-2/test-2-rnd-polygons-1-4.off", ks, num_iters, results, all_times, num_tests));
  results = {7,1,14,24,13,2};
  assert(run_test<Traits>("data/stress-test-2/test-3-rnd-polygons-1-4.off", ks, num_iters, results, all_times, num_tests));
  results = {7,1,14,24,13,2};
  assert(run_test<Traits>("data/stress-test-2/test-4-rnd-polygons-1-3.off", ks, num_iters, results, all_times, num_tests));
  results = {8,1,19,35,20,3};
  assert(run_test<Traits>("data/stress-test-2/test-5-rnd-polygons-2-4.off", ks, num_iters, results, all_times, num_tests));
  results = {9,1,26,50,30,5};
  assert(run_test<Traits>("data/stress-test-2/test-6-rnd-polygons-3-4.off", ks, num_iters, results, all_times, num_tests));*/

  // Stress tests 3.
/*
  results = {8,1,20,37,21,3};
  assert(run_test<Traits>("data/stress-test-3/test-1-rnd-polygons-2-3.off" , ks, num_iters, results, all_times, num_tests));
  results = {8,1,17,30,17,3};
  assert(run_test<Traits>("data/stress-test-3/test-2-rnd-polygons-2-3.off" , ks, num_iters, results, all_times, num_tests));
  results = {8,1,19,35,20,3};
  assert(run_test<Traits>("data/stress-test-3/test-3-rnd-polygons-2-3.off" , ks, num_iters, results, all_times, num_tests));
  results = {8,1,19,35,20,3};
  assert(run_test<Traits>("data/stress-test-3/test-4-rnd-polygons-2-4.off" , ks, num_iters, results, all_times, num_tests));
  results = {7,1,10,18,11,2};
  assert(run_test<Traits>("data/stress-test-3/test-5-rnd-polygons-1-3.off" , ks, num_iters, results, all_times, num_tests));
  results = {8,1,19,35,20,3};
  assert(run_test<Traits>("data/stress-test-3/test-6-rnd-polygons-2-3.off" , ks, num_iters, results, all_times, num_tests));
  results = {8,1,22,41,23,3};
  assert(run_test<Traits>("data/stress-test-3/test-7-rnd-polygons-2-4.off" , ks, num_iters, results, all_times, num_tests));
  results = {8,1,18,33,19,3};
  assert(run_test<Traits>("data/stress-test-3/test-8-rnd-polygons-2-10.off", ks, num_iters, results, all_times, num_tests));
  results = {10,1,39,82,50,7};
  assert(run_test<Traits>("data/stress-test-3/test-9-rnd-polygons-4-4.off" , ks, num_iters, results, all_times, num_tests));
  results = {11,1,55,119,78,13};
  assert(run_test<Traits>("data/stress-test-3/test-10-rnd-polygons-5-4.off", ks, num_iters, results, all_times, num_tests));*/

  // Stress tests 4.
/*
  results = {8,1,20,37,21,3};
  assert(run_test<Traits>("data/stress-test-4/test-1-rnd-polygons-2-6.off" , ks, num_iters, results, all_times, num_tests));
  results = {9,1,29,58,36,6};
  assert(run_test<Traits>("data/stress-test-4/test-2-rnd-polygons-3-8.off" , ks, num_iters, results, all_times, num_tests));
  results = {10,1,37,76,48,8};
  assert(run_test<Traits>("data/stress-test-4/test-3-rnd-polygons-4-4.off" , ks, num_iters, results, all_times, num_tests));
  results = {10,1,37,77,46,6};
  assert(run_test<Traits>("data/stress-test-4/test-4-rnd-polygons-4-6.off" , ks, num_iters, results, all_times, num_tests));
  results = {12,2,83,191,129,21};
  assert(run_test<Traits>("data/stress-test-4/test-5-rnd-polygons-6-4.off" , ks, num_iters, results, all_times, num_tests));
  results = {11,1,50,107,71,14};
  assert(run_test<Traits>("data/stress-test-4/test-6-rnd-polygons-5-6.off" , ks, num_iters, results, all_times, num_tests));
  results = {13,2,104,246,160,23};
  assert(run_test<Traits>("data/stress-test-4/test-7-rnd-polygons-7-6.off" , ks, num_iters, results, all_times, num_tests));
  results = {13,1,69,152,96,13};
  assert(run_test<Traits>("data/stress-test-4/test-8-rnd-polygons-7-8.off" , ks, num_iters, results, all_times, num_tests));
  results = {18,3,250,629,449,76};
  assert(run_test<Traits>("data/stress-test-4/test-9-rnd-polygons-12-4.off", ks, num_iters, results, all_times, num_tests));*/

  // Stress tests 5.

  results = {21,2,468,1224,720,66};
  run_test<Traits>("data/stress-test-5/test-1-rnd-polygons-15-6.off", ks, num_iters, results, all_times, num_tests);
  //results = {26,3,1037,2829,1693,161};
  //run_test<Traits>("data/stress-test-5/test-2-rnd-polygons-20-4.off", ks, num_iters, results, all_times, num_tests);

  // Real data tests.
/*
  results = {16,1,133,315,212,34};
  assert(run_test<Traits>("data/real-data-test/test-10-polygons.off", ks, num_iters, results, all_times, num_tests));
  results = {18,2,217,543,370,58};
  run_test<Traits>("data/real-data-test/test-15-polygons.off", ks, num_iters, results, all_times, num_tests);
  //results = {21,3,375,974,629,74};
  //run_test<Traits>("data/real-data-test/test-20-polygons.off", ks, num_iters, results, all_times, num_tests);*/

  //ks.clear();
  //ks.push_back(5);

  //results = { 38,3,2556,7128,3272,133 }; // fails for k = 1 and coplanarity = 0.1; and k = 6 and coplanarity = 0.5
  //run_test<Traits>("data/real-data-test/test-40-polygons.ply", ks, num_iters, results, all_times, num_tests);

  std::cout << std::endl << "--OUTPUT STATS:" << std::endl;
  std::cout << "* number of tests: "               << num_tests << std::endl;
  std::cout << "* number of iterations per test: " << num_iters << std::endl;
  std::cout << "* k intersections: {";
  for (const auto k : ks) {
    std::cout << k << ",";
  }
  std::cout << "...}" << std::endl;

  if (num_tests != 0) {
    std::cout << std::endl << "--TIMINGS:" << std::endl;
    for (std::size_t i = 0; i < all_times.size(); ++i) {
      std::cout << "* time (sec.), test #" << std::to_string(i) << ": {";
      for (const double& time : all_times[i]) {
        std::cout << time << ", ";
      }
      std::cout << "...}" << std::endl;
    }
  }

  const auto kernel_name = boost::typeindex::type_id<typename Traits::Kernel>().pretty_name();
  const auto intersection_kernel_name = boost::typeindex::type_id<typename Traits::Intersection_Kernel>().pretty_name();
  if (num_tests != 0) {
    std::cout << std::endl << kernel_name << " with " << intersection_kernel_name << " intersections" <<
    ": ALL " << num_tests << " TESTS SUCCESS!" << std::endl << std::endl;
  }
  else {
    std::cout << std::endl << kernel_name << " with " << intersection_kernel_name << " intersections" <<
    ": ALL " << num_tests << " TESTS FAILED!" << std::endl << std::endl;
  }
}

#include <CGAL/intersections.h>

int main(const int /* argc */, const char** /* argv */) {
  // run_all_tests<SCF>();
  // run_all_tests<SCD>();
  //run_all_tests<EPECK>();

  // Passes all tests except for those when
  // intersections lead to accumulated errors.
  run_all_tests<Traits>();
  return EXIT_SUCCESS;
}
