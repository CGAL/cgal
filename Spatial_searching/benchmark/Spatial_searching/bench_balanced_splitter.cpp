#include <cmath>
#include <vector>
#include <fstream>

#include <CGAL/Real_timer.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/centroid.h>
#include <CGAL/Kd_tree.h>

using SCF   = CGAL::Simple_cartesian<float>;
using SCD   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using Timer = CGAL::Real_timer;

template<typename Kernel>
void save_bbox(
  const std::string name,
  const std::array<typename Kernel::Point_3, 8>& bbox) {

  using Point_3 = typename Kernel::Point_3;
  using Surface_mesh = CGAL::Surface_mesh<Point_3>;
  Surface_mesh surface_mesh;
  CGAL::make_hexahedron(
    bbox[0], bbox[1], bbox[2], bbox[3],
    bbox[4], bbox[5], bbox[6], bbox[7],
    surface_mesh);
  CGAL::IO::write_polygon_mesh(name, surface_mesh);
}

template<typename Kernel>
bool is_centered(std::vector<typename Kernel::Point_3>& point_set) {

  using FT = typename Kernel::FT;
  const auto centroid = CGAL::centroid(point_set.begin(), point_set.end());
  const FT tol  = FT(1) / FT(1000000);
  return (
    CGAL::abs(centroid.x()) < tol &&
    CGAL::abs(centroid.y()) < tol &&
    CGAL::abs(centroid.z()) < tol );
}

template<typename Kernel>
typename Kernel::Point_3 translate_point_set(
  const typename Kernel::Vector_3& translation,
  std::vector<typename Kernel::Point_3>& point_set) {

  using Affine_transformation_3 = CGAL::Aff_transformation_3<Kernel>;
  for (auto& point : point_set) {
    point = point.transform(
      Affine_transformation_3(CGAL::Translation(), translation));
  }
  std::ofstream out("1-translated.xyz");
  CGAL::IO::write_XYZ(out, point_set);
  return CGAL::centroid(point_set.begin(), point_set.end());
}

template<typename Kernel>
void initialize_bbox(
  const std::vector<typename Kernel::Point_3>& point_set,
  std::array<typename Kernel::Point_3, 8>& bbox) {

  using Point_3 = typename Kernel::Point_3;
  const auto box = CGAL::bbox_3(point_set.begin(), point_set.end());

  // The order of faces corresponds to the standard order from here:
  // https://doc.cgal.org/latest/BGL/group__PkgBGLHelperFct.html#gad9df350e98780f0c213046d8a257358e
  bbox = {
    Point_3(box.xmin(), box.ymin(), box.zmin()),
    Point_3(box.xmax(), box.ymin(), box.zmin()),
    Point_3(box.xmax(), box.ymax(), box.zmin()),
    Point_3(box.xmin(), box.ymax(), box.zmin()),
    Point_3(box.xmin(), box.ymax(), box.zmax()),
    Point_3(box.xmin(), box.ymin(), box.zmax()),
    Point_3(box.xmax(), box.ymin(), box.zmax()),
    Point_3(box.xmax(), box.ymax(), box.zmax()) };

  save_bbox<Kernel>("2-tight-bbox.off", bbox);
}

template<typename Kernel>
typename Kernel::FT enlarge_bounding_box(
  const typename Kernel::FT enlarge_bbox_ratio,
  std::array<typename Kernel::Point_3, 8>& bbox) {

  using Transform_3 = typename Kernel::Aff_transformation_3;
  const auto a = CGAL::centroid(bbox.begin(), bbox.end());
  Transform_3 scale(CGAL::Scaling(), enlarge_bbox_ratio);
  for (auto& point : bbox)
    point = scale.transform(point);

  const auto b = CGAL::centroid(bbox.begin(), bbox.end());
  Transform_3 translate(CGAL::Translation(), a - b);
  for (auto& point : bbox)
    point = translate.transform(point);

  const auto& minp = bbox[0];
  const auto& maxp = bbox[7];
  const auto side_1 = maxp.x() - minp.x();
  const auto side_2 = maxp.y() - minp.y();
  const auto side_3 = maxp.z() - minp.z();

  save_bbox<Kernel>("3-enlarged-bbox.off", bbox);
  return (CGAL::max)((CGAL::max)(side_1, side_2), side_3);
}

template<typename Kernel>
void generate_query_points(
  const typename Kernel::FT a,
  const std::size_t num_query_points,
  const std::array<typename Kernel::Point_3, 8>& bbox,
  std::vector<typename Kernel::Point_3>& query_points) {

  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;

  query_points.clear();
  query_points.reserve(num_query_points);
  const std::size_t total_size = num_query_points * 2;
  CGAL::Random_points_in_cube_3<Point_3> generator(a / FT(2));

  const auto& minp = bbox[0];
  const auto& maxp = bbox[7];

  std::size_t count = 0;
  for (std::size_t i = 0; i < total_size; ++i) {
    const auto query_point = *generator++;
    if (
      query_point.x() >= minp.x() && query_point.x() <= maxp.x() &&
      query_point.y() >= minp.y() && query_point.y() <= maxp.y() &&
      query_point.z() >= minp.z() && query_point.z() <= maxp.z() ) {

      query_points.push_back(query_point);
      ++count;
    }
    if (count == num_query_points) break;
  }
  assert(query_points.size() == num_query_points);

  std::ofstream out("4-query-points.xyz");
  CGAL::IO::write_XYZ(out, query_points);
}

template<typename Kernel, typename Traits, typename Splitter>
void bench_splitter(
  const std::size_t bucket_size,
  const std::size_t k,
  const std::vector<typename Kernel::Point_3>& point_set,
  const std::vector<typename Kernel::Point_3>& query_points) {

  std::cout << "- num input points: " << point_set.size() << std::endl;
  std::cout << "- num query points: " << query_points.size() << std::endl;
  std::cout << "- bucket size: " << bucket_size << std::endl;
  std::cout << "- num k: " << k << std::endl;

  Timer timer;
  using Kd_tree = CGAL::Kd_tree<Traits, Splitter>;
  using Distance = CGAL::Euclidean_distance<Traits>;
  using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Traits, Distance, Splitter, Kd_tree>;

  // Building.
  Splitter splitter(bucket_size);
  Kd_tree tree(point_set.begin(), point_set.end(), splitter);
  std::cout << "- building the tree ... ";
  timer.reset();
  timer.start();
  tree.build();
  timer.stop();
  std::cout << "done in: " << timer.time() << " sec." << std::endl;

  // K neighbor search.
  std::cout << "- applying k neighbor search ... ";
  timer.reset();
  timer.start();
  for (const auto& query_point : query_points) {
    Neighbor_search nsearch(tree, query_point, k);
    long count = 0;
    for (auto it = nsearch.begin(); it != nsearch.end(); ++it) {
      ++count;
    }
    assert(std::distance(nsearch.begin(), nsearch.end()) == count);
  }
  timer.stop();
  std::cout << "done in: " << timer.time() << " sec." << std::endl;
}

template<typename Kernel>
void bench_random_in_bbox(
  const std::string filename,
  const std::size_t num_query_points = 100,
  const std::size_t bucket_size = 10,
  const std::size_t k = 6) {

  std::cout.precision(20);
  std::cout << std::endl;
  std::cout << "---- RANDOM IN BBOX BENCH" << std::endl;
  std::cout << "- filename " << filename << " ... " << std::endl;

  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Vector_3 = typename Kernel::Vector_3;
  using Traits = CGAL::Search_traits_3<Kernel>;
  using Balanced_splitter = CGAL::Balanced_splitter<Traits>;
  using Sliding_midpoint = CGAL::Sliding_midpoint<Traits>;

  // Read point set.
  std::cout << std::endl << "* reading input point set" << std::endl;
  std::vector<Point_3> point_set;
  std::ifstream in(filename);
  CGAL::IO::set_ascii_mode(in);
  Point_3 p;
  while (in >> p) {
    point_set.push_back(p);
  }
  std::cout << "- number of points: " << point_set.size() << std::endl;
  auto centroid = CGAL::centroid(point_set.begin(), point_set.end());
  std::cout << "- centroid: " << centroid << std::endl;
  std::ofstream out("0-original.xyz");
  CGAL::IO::write_XYZ(out, point_set);
  assert(point_set.size() > 0);

  // Translate point set.
  std::cout << std::endl << "* translating input point set" << std::endl;
  const Vector_3 translation(centroid, CGAL::ORIGIN);
  centroid = translate_point_set<Kernel>(translation, point_set);
  std::cout << "- centroid: " << centroid << std::endl;
  assert(is_centered<Kernel>(point_set));

  // Create a bbox.
  std::cout << std::endl << "* creating tight bbox" << std::endl;
  std::array<Point_3, 8> bbox;
  initialize_bbox<Kernel>(point_set, bbox);
  std::cout << "- axis aligned bbox with 8 vertices" << std::endl;

  // Enlarge the bbox.
  std::cout << std::endl << "* enlarging bbox" << std::endl;
  const FT enlarge_bbox_ratio = FT(3) / FT(2);
  std::cout << "- enlarge bbox ratio: " << enlarge_bbox_ratio << std::endl;
  const FT longest_side = enlarge_bounding_box<Kernel>(enlarge_bbox_ratio, bbox);
  std::cout << "- longest side: " << longest_side << std::endl;

  // Generating random points.
  std::cout << std::endl << "* generating random points" << std::endl;
  std::cout << "- num queries: " << num_query_points << std::endl;
  std::vector<Point_3> query_points;
  generate_query_points<Kernel>(longest_side, num_query_points, bbox, query_points);
  std::cout << "- num random points in bbox: " << query_points.size() << std::endl;

  // Benching the tree.
  std::cout << std::endl << "* processing balanced splitter" << std::endl;
  bench_splitter<Kernel, Traits, Balanced_splitter>(bucket_size, k, point_set, query_points);

  std::cout << std::endl << "* processing sliding midpoint splitter" << std::endl;
  bench_splitter<Kernel, Traits, Sliding_midpoint>(bucket_size, k, point_set, query_points);

  std::cout << std::endl;
}

template<typename Kernel>
void bench_translated(
  const std::string filename,
  const typename Kernel::FT d = 0.1,
  const std::size_t bucket_size = 10,
  const std::size_t k = 6) {

  std::cout.precision(20);
  std::cout << std::endl;
  std::cout << "---- TRANSLATION BENCH" << std::endl;
  std::cout << "- filename " << filename << " ... " << std::endl;

  using Point_3 = typename Kernel::Point_3;
  using Vector_3 = typename Kernel::Vector_3;
  using Traits = CGAL::Search_traits_3<Kernel>;
  using Balanced_splitter = CGAL::Balanced_splitter<Traits>;
  using Sliding_midpoint = CGAL::Sliding_midpoint<Traits>;

  // Read point set.
  std::cout << std::endl << "* reading input point set" << std::endl;
  std::vector<Point_3> point_set;
  std::ifstream in(filename);
  CGAL::IO::set_ascii_mode(in);
  Point_3 p;
  while (in >> p) {
    point_set.push_back(p);
  }
  std::cout << "- number of points: " << point_set.size() << std::endl;
  auto centroid = CGAL::centroid(point_set.begin(), point_set.end());
  std::cout << "- centroid: " << centroid << std::endl;
  std::ofstream out("0-original.xyz");
  CGAL::IO::write_XYZ(out, point_set);
  assert(point_set.size() > 0);

  // Translate point set.
  std::cout << std::endl << "* creating query points by translation" << std::endl;
  auto query_points = point_set;
  const Vector_3 translation(d, d, d);
  centroid = translate_point_set<Kernel>(translation, query_points);
  std::cout << "- centroid: " << centroid << std::endl;
  std::cout << "- num query points: " << query_points.size() << std::endl;

  // Benching the tree.
  std::cout << std::endl << "* processing balanced splitter" << std::endl;
  std::cout << "- translation distance: " << d << std::endl;
  bench_splitter<Kernel, Traits, Balanced_splitter>(bucket_size, k, point_set, query_points);

  std::cout << std::endl << "* processing sliding midpoint splitter" << std::endl;
  std::cout << "- translation distance: " << d << std::endl;
  bench_splitter<Kernel, Traits, Sliding_midpoint>(bucket_size, k, point_set, query_points);
}

int main(int argc, char* argv[]) {

  double d;
  std::size_t num_query_points, bucket_size, k;
  const std::string filename = (argc > 1 ? argv[1] : "data/sphere.xyz");

  // BENCH 1.
  num_query_points = 100;
  bucket_size = 10;
  k = 6;
  bench_random_in_bbox<SCD>(filename, num_query_points, bucket_size, k);

  // BENCH 2.
  d = 0.05;
  bucket_size = 10;
  k = 1;
  bench_translated<SCD>(filename, d, bucket_size, k);
}
