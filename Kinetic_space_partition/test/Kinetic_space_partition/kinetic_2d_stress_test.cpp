#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>

#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Kinetic_space_partition_2.h>

using SC = CGAL::Simple_cartesian<double>;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK_to_EPICK = CGAL::Cartesian_converter<EPECK, EPICK>;

using FT        = typename EPECK::FT;
using Point_2   = typename EPECK::Point_2;
using Vector_2  = typename EPECK::Vector_2;
using Segment_2 = typename EPECK::Segment_2;
using Transform = CGAL::Aff_transformation_2<EPECK>;

using Random = CGAL::Random;
using Mesh   = CGAL::Surface_mesh<Point_2>;

using Exact_reconstruction   = CGAL::Kinetic_space_partition_2<EPECK>;
using Inexact_reconstruction = CGAL::Kinetic_space_partition_2<EPICK>;

Random cgal_rand;

void add_regular_case(
  std::vector<Segment_2>& segments) {

  const std::size_t size_before = segments.size();
  segments.push_back(Segment_2(Point_2(0.0, 1), Point_2(0.0, 3)));
  segments.push_back(Segment_2(Point_2(0.0, 5), Point_2(0.0, 7)));
  segments.push_back(Segment_2(Point_2(4.0, 1), Point_2(4.0, 3)));
  segments.push_back(Segment_2(Point_2(4.0, 6), Point_2(4.0, 7)));
  segments.push_back(Segment_2(Point_2(1.0, 0), Point_2(3.0, 0)));
  segments.push_back(Segment_2(Point_2(2.0, 4), Point_2(3.0, 4)));
  segments.push_back(Segment_2(Point_2(1.2, 8), Point_2(2.5, 8)));

  // Random rotation.
  const double sine = cgal_rand.get_double(-1, 1);
  const double cosine = CGAL::sqrt(1.0 - sine * sine);

  Transform rotate(CGAL::Rotation(), sine, cosine);
  Transform scale(CGAL::Scaling(), cgal_rand.get_double(0.1, 10));
  Transform translate(CGAL::Translation(), Vector_2(
    cgal_rand.get_double(-5, 5), cgal_rand.get_double(-5, 5)));
  Transform transform = scale * rotate * translate;

  for (std::size_t i = size_before; i < segments.size(); ++i) {
    const auto source = transform.transform(segments[i].source());
    const auto target = transform.transform(segments[i].target());
    segments[i] = Segment_2(source, target);
  }
}

void add_star_case(
  std::vector<Segment_2>& segments,
  std::size_t star_branches) {

  const std::size_t size_before = segments.size();
  Segment_2 base(Point_2(0, 1), Point_2(0, 3));

  for (std::size_t i = 0; i < star_branches; ++i) {
    const double angle = 2.0 * CGAL_PI * (i / double(star_branches));
    Transform rotate(CGAL::Rotation(), std::sin(angle), std::cos(angle));
    segments.push_back(Segment_2(
      rotate.transform(base.source()),
      rotate.transform(base.target())));
  }

  // Random rotation.
  const double sine = cgal_rand.get_double(-1.1);
  const double cosine = CGAL::sqrt(1.0 - sine * sine);
  Transform rotate(CGAL::Rotation(), sine, cosine);
  Transform scale(CGAL::Scaling(), cgal_rand.get_double(0.1, 10));
  Transform translate(CGAL::Translation(), Vector_2(
    cgal_rand.get_double(-5, 5), cgal_rand.get_double(-5, 5)));
  Transform transform = scale * rotate * translate;

  for (std::size_t i = size_before; i < segments.size(); ++i) {
    const auto source = transform.transform(segments[i].source());
    const auto target = transform.transform(segments[i].target());
    segments[i] = Segment_2(source, target);
  }
}

template<typename Kernel>
void get_segments_from_exact(
  const std::vector<Segment_2>&,
  std::vector<typename Kernel::Segment_2>&) {
  abort();
}

template<>
void get_segments_from_exact<EPECK>(
  const std::vector<Segment_2>& exact_segments,
  std::vector<Segment_2>& segments) {

  segments.reserve(exact_segments.size());
  std::copy(exact_segments.begin(), exact_segments.end(), std::back_inserter(segments));
}

template<>
void get_segments_from_exact<EPICK>(
  const std::vector<Segment_2>& exact_segments,
  std::vector<typename EPICK::Segment_2>& segments) {

  static EPECK_to_EPICK e2e;
  segments.reserve(exact_segments.size());
  std::transform(exact_segments.begin(), exact_segments.end(),
    std::back_inserter(segments),
    [&](const Segment_2& segment) -> typename EPICK::Segment_2 {
      return e2e(segment);
    });
}

template<typename Kernel>
void test_segments(
  std::string test_name,
  const std::vector<Segment_2>& exact_segments,
  unsigned int k) {

  CGAL::Real_timer t;
  t.start();

  std::vector<typename Kernel::Segment_2> segments;
  get_segments_from_exact<Kernel>(exact_segments, segments);

  CGAL::Kinetic_space_partition_2<Kernel> ksp;
  ksp.partition(segments, CGAL::Identity_property_map<typename Kernel::Segment_2>(), k, 2);
  segments.clear();
  ksp.output_partition_edges_to_segment_soup(std::back_inserter(segments));

#ifdef OUTPUT_FILES
  std::ofstream output_file(
    test_name + (std::is_same<Kernel, EPECK>::value ? "_exact" : "_inexact") + "_output.polylines.txt");
  for (const auto& s : segments) {
    output_file << "2 " << s.source() << " 0 " << s.target() << " 0" << std::endl;
  }
#endif

  if (!ksp.check_integrity(true)) {
    std::cerr << "ERROR: Integrity of reconstruction failed!" << std::endl;
    return;
  }

  CGAL::Surface_mesh<typename Kernel::Point_2> mesh;
  if (ksp.output_partition_cells_to_face_graph(mesh)) {
#ifdef OUTPUT_FILES
    std::ofstream output_shapes_file(
      test_name + (std::is_same<Kernel, EPECK>::value ? "_exact" : "_inexact") + "_faces.ply");
    output_shapes_file
      << "ply" << std::endl
      << "format ascii 1.0" << std::endl
      << "element vertex " << mesh.number_of_vertices() << std::endl
      << "property double x" << std::endl
      << "property double y" << std::endl
      << "property double z" << std::endl
      << "element face " << mesh.number_of_faces() << std::endl
      << "property list uchar int vertex_index" << std::endl
      << "property uchar red" << std::endl
      << "property uchar green" << std::endl
      << "property uchar blue" << std::endl
      << "end_header" << std::endl;

    for (const auto& vindex : vertices(mesh)) {
      output_shapes_file << mesh.point(vindex) << " 0" << std::endl;
    }
    for (const auto& findex : faces(mesh)) {
      output_shapes_file << degree(findex, mesh);
      for (const auto& hindex : CGAL::halfedges_around_face(halfedge(findex, mesh), mesh)) {
        output_shapes_file << " " << int(target(hindex,mesh));
      }
      output_shapes_file
        << " " << cgal_rand.get_int(64,192)
        << " " << cgal_rand.get_int(64,192)
        << " " << cgal_rand.get_int(64,192) << std::endl;
    }
#endif
  }
  else {
    std::cerr << "ERROR: Invalid face graph!" << std::endl;
  }

  t.stop();
  std::cerr
    << " -> "
    << (std::is_same<Kernel, EPECK>::value ? "exact " : "inexact ")
    << "stress test " << test_name << " done in " << t.time() << " seconds" << std::endl;
}

void stress_test(
  std::string test_name,
  std::size_t nb_random_lines,
  std::size_t nb_regular_boxes,
  std::size_t nb_stars,
  std::size_t star_branches,
  std::size_t k) {

  cgal_rand = CGAL::Random(0);

  std::cerr << "[Stress test " << test_name << "]" << std::endl;
  std::vector<Segment_2> exact_segments;

  for (std::size_t i = 0; i < nb_regular_boxes; ++i) {
    add_regular_case(exact_segments);
  }

  for (std::size_t i = 0; i < nb_stars; ++i) {
    add_star_case(exact_segments, star_branches);
  }

  CGAL::Bbox_2 bbox(0, 0, 5, 5);
  if (!exact_segments.empty()) {
    for (const Segment_2& segment : exact_segments) {
      bbox = bbox + segment.bbox();
    }
  }

  Point_2 pmin(bbox.xmin(), bbox.ymin());
  Point_2 pmax(bbox.xmax(), bbox.ymax());
  double seg_size = CGAL::to_double(FT(0.1) *
    CGAL::approximate_sqrt(CGAL::squared_distance(pmin, pmax)));

  for (std::size_t i = 0; i < nb_random_lines; ++i) {

    Point_2 source(
      cgal_rand.get_double(bbox.xmin(), bbox.xmax()),
      cgal_rand.get_double(bbox.ymin(), bbox.ymax()));
    Vector_2 vec(
      cgal_rand.get_double(-seg_size, seg_size),
      cgal_rand.get_double(-seg_size, seg_size));
    Point_2 target = source + vec;
    exact_segments.push_back(Segment_2(source, target));
  }

#ifdef OUTPUT_FILES
  std::ofstream input_file(test_name + "_input.polylines.txt");
  for (const Segment_2& s : exact_segments) {
    input_file << "2 " << s.source() << " 0 " << s.target() << " 0" << std::endl;
  }
#endif

#ifdef TEST_EPECK
  if (exact_segments.size() < 500) {
    test_segments<EPECK>(test_name, exact_segments, k);
  }
  else {
    std::cerr << " -> skipping exact test to avoid overly long running time (too many segments)" << std::endl;
  }
#endif

#if true
  test_segments<EPICK>(test_name, exact_segments, k);
#endif
}

int main(const int /* argc */, const char** /* argv */) {

  CGAL::Real_timer t;
  t.start();

  stress_test("01_30_random_lines", 30, 0, 0, 0, 2);
  stress_test("02_300_random_lines", 300, 0, 0, 0, 2);
  stress_test("03_300_random_lines_k_10", 300, 0, 0, 0, 10);
#if true
  stress_test("04_3000_random_lines", 3000, 0, 0, 0, 2);
  stress_test("05_3000_random_lines_k_3", 3000, 0, 0, 0, 3);
#endif

  stress_test("06_regular_case", 0, 1, 0, 0, 2);
  stress_test("07_multi_regular_case", 0, 5, 0, 0, 2);
  stress_test("08_multi_regular_case_and_random_lines", 30, 5, 0, 0, 2);
  stress_test("09_big_multi_regular_case_and_random_lines", 100, 30, 0, 0, 4);

  stress_test("10_cross", 0, 0, 1, 4, 2);
  stress_test("11_star", 0, 0, 1, 6, 2);
  stress_test("12_multiple_stars", 0, 0, 5, 12, 2);
  stress_test("13_stars_and_regular", 0, 5, 5, 12, 3);
  stress_test("14_everything", 100, 30, 5, 12, 2);
#if true
  stress_test("15_mayhem", 3000, 100, 10, 20, 4);
#endif

  t.stop();
  std::cerr << "All tests done in " << t.time() << " seconds!" << std::endl;
  return EXIT_SUCCESS;
}
