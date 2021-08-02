#include <cmath>
#include <vector>
#include <fstream>

#include <boost/type_index.hpp>

#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/centroid.h>
#include <CGAL/Kd_tree.h>

using Timer = CGAL::Real_timer;

using SCF   = CGAL::Simple_cartesian<float>;
using SCD   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

using SCF_traits   = CGAL::Search_traits_3<SCF>;
using SCD_traits   = CGAL::Search_traits_3<SCD>;
using EPICK_traits = CGAL::Search_traits_3<EPICK>;
using EPECK_traits = CGAL::Search_traits_3<EPECK>;

using Color = std::array<unsigned char, 3>;

// Define how a color should be stored.
namespace CGAL {
  template<class F>
  struct Output_rep< ::Color, F > {
    const ::Color& c;
    static const bool is_specialized = true;
    Output_rep(const ::Color& c) : c(c) { }
    std::ostream& operator()(std::ostream& out) const {
      if (IO::is_ascii(out)) {
        out << int(c[0]) << " " << int(c[1]) << " " << int(c[2]);
      } else {
        out.write(reinterpret_cast<const char*>(&c), sizeof(c));
      }
      return out;
    }
  };
} // namespace CGAL

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
  const bool verbose,
  const typename Kernel::Vector_3& translation,
  std::vector<typename Kernel::Point_3>& point_set) {

  using Affine_transformation_3 = CGAL::Aff_transformation_3<Kernel>;
  for (auto& point : point_set) {
    point = point.transform(
      Affine_transformation_3(CGAL::Translation(), translation));
  }

  if (verbose) {
    std::ofstream outfile("1-translated.xyz");
    CGAL::IO::write_XYZ(outfile, point_set);
  }
  return CGAL::centroid(point_set.begin(), point_set.end());
}

template<typename Kernel>
void initialize_bounding_box(
  const bool verbose,
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

  if (verbose) {
    save_bbox<Kernel>("2-tight-bbox.off", bbox);
  }
}

template<typename Kernel>
typename Kernel::FT enlarge_bounding_box(
  const bool verbose,
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

  if (verbose) {
    save_bbox<Kernel>("3-enlarged-bbox.off", bbox);
  }
  return (CGAL::max)((CGAL::max)(side_1, side_2), side_3);
}

template<typename Kernel>
void generate_query_points(
  const bool verbose,
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

  if (verbose) {
    std::ofstream outfile("4-query-points.xyz");
    CGAL::IO::write_XYZ(outfile, query_points);
  }
}

template<typename Kernel, typename Traits, typename Splitter>
std::pair<double, double> bench_splitter(
  const bool verbose,
  const std::size_t num_iters,
  const std::size_t bucket_size,
  const std::size_t k,
  const std::vector<typename Kernel::Point_3>& point_set,
  const std::vector<typename Kernel::Point_3>& query_points,
  const bool save_distribution = false) {

  Timer timer;
  using Point_3 = typename Kernel::Point_3;
  using Kd_tree = CGAL::Kd_tree<Traits, Splitter>;
  using Distance = CGAL::Euclidean_distance<Traits>;
  using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Traits, Distance, Splitter, Kd_tree>;

  // Building the tree.
  double avg_build_time = 0.0;
  Splitter splitter(bucket_size);
  Kd_tree tree(splitter);
  for (std::size_t iter = 0; iter < num_iters; ++iter) {
    tree.clear();
    tree.insert(point_set.begin(), point_set.end());
    timer.reset();
    timer.start();
    tree.build();
    timer.stop();
    avg_build_time += timer.time();
  }
  CGAL_assertion(num_iters > 0);
  avg_build_time /= static_cast<double>(num_iters);

  if (verbose) {
    std::cout << "- AVG BUILDS TIME: " << avg_build_time << " sec." << std::endl;

    if (save_distribution) {
      using Point_with_color = std::pair<Point_3, Color>;
      using PLY_Point_map = CGAL::First_of_pair_property_map<Point_with_color>;
      using PLY_Color_map = CGAL::Second_of_pair_property_map<Point_with_color>;
      using Point_with_index = std::pair<Point_3, std::size_t>;

      std::ofstream outfile("5-leaves-distribution.ply");
      std::vector<Point_with_index> pwi;
      tree.print(std::back_inserter(pwi));
      CGAL_assertion(pwi.size() > 0);

      std::vector<Point_with_color> pwc;
      pwc.reserve(pwi.size());
      for (const auto& pair : pwi) {
        CGAL::Random rand(static_cast<unsigned int>(pair.second));
        const unsigned char r =
          static_cast<unsigned char>(64 + rand.get_int(0, 192));
        const unsigned char g =
          static_cast<unsigned char>(64 + rand.get_int(0, 192));
        const unsigned char b =
          static_cast<unsigned char>(64 + rand.get_int(0, 192));
        const Color color = CGAL::make_array(r, g, b);
        pwc.push_back(std::make_pair(pair.first, color));
      }
      CGAL_assertion(pwc.size() == pwi.size());

      CGAL::IO::set_ascii_mode(outfile);
      CGAL::IO::write_PLY_with_properties(
        outfile, pwc,
        CGAL::make_ply_point_writer(PLY_Point_map()),
          std::make_tuple(
            PLY_Color_map(),
            CGAL::PLY_property<unsigned char>("red"),
            CGAL::PLY_property<unsigned char>("green"),
            CGAL::PLY_property<unsigned char>("blue")));
    }
  }

  // K neighbor search.
  double avg_search_time = 0.0;
  for (std::size_t iter = 0; iter < num_iters; ++iter) {
    timer.reset();
    timer.start();
    for (const auto& query_point : query_points) {
      CGAL_assertion(tree.size() != 0);
      Neighbor_search nsearch(tree, query_point, k);
      long count = 0;
      for (auto it = nsearch.begin(); it != nsearch.end(); ++it) {
        ++count;
      }
      assert(std::distance(nsearch.begin(), nsearch.end()) == count);
    }
    timer.stop();
    avg_search_time += timer.time();
  }
  CGAL_assertion(num_iters > 0);
  avg_search_time /= static_cast<double>(num_iters);

  if (verbose) {
    std::cout << "- AVG SEARCH TIME: " << avg_search_time << " sec." << std::endl;
  }
  return std::make_pair(avg_build_time, avg_search_time);
}

template<typename Kernel, typename Traits, typename Splitter>
std::vector<double> bench_random_in_bbox(
  const std::string filename, const bool verbose,
  std::size_t num_query_points, std::size_t num_iters,
  std::size_t bucket_size, std::size_t k) {

  if (verbose) {
    std::cout << std::endl;
    std::cout << "---- RANDOM IN BBOX BENCH" << std::endl;
    std::cout << "- filename " << filename << " ... " << std::endl << std::endl;
  }

  std::vector<double> out;
  out.push_back(k);
  out.push_back(bucket_size);

  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Vector_3 = typename Kernel::Vector_3;

  // Read point set.
  Timer timer;
  timer.start();
  std::vector<Point_3> point_set;
  std::ifstream in(filename);
  CGAL::IO::set_ascii_mode(in);
  Point_3 p;
  while (in >> p) {
    point_set.push_back(p);
  }
  auto centroid = CGAL::centroid(point_set.begin(), point_set.end());
  if (num_query_points == 0) num_query_points = point_set.size();
  assert(point_set.size() > 0);
  timer.stop();
  out.push_back(timer.time());

  if (verbose) {
    std::cout << "- PARAMETERS: " << std::endl;
    std::cout << "-            num k: " << k << std::endl;
    std::cout << "-      bucket size: " << bucket_size << std::endl;
    std::cout << "- num input points: " << point_set.size() << std::endl;
    std::cout << "- num query points: " << num_query_points << std::endl;
    std::cout << std::endl;

    std::cout << "- CREATE POINT TIME: " << out.back() << " sec." << std::endl;

    std::ofstream outfile("0-original.xyz");
    CGAL::IO::write_XYZ(outfile, point_set);
  }

  // Translate point set.
  const Vector_3 translation(centroid, CGAL::ORIGIN);
  translate_point_set<Kernel>(verbose, translation, point_set);
  assert(is_centered<Kernel>(point_set));

  // Create tight bbox.
  std::array<Point_3, 8> bbox;
  initialize_bounding_box<Kernel>(verbose, point_set, bbox);

  // Enlarge bbox.
  const FT enlarge_bbox_ratio = FT(3) / FT(2);
  const FT longest_side = enlarge_bounding_box<Kernel>(
    verbose, enlarge_bbox_ratio, bbox);

  // Generate random points in bbox.
  timer.reset();
  timer.start();
  std::vector<Point_3> query_points;
  generate_query_points<Kernel>(
    verbose, longest_side, num_query_points, bbox, query_points);
  timer.stop();
  out.push_back(timer.time());

  if (verbose) {
    std::cout << "- CREATE QUERY TIME: " << out.back() << " sec." << std::endl;
  }

  // Bench tree.
  if (verbose) std::cout << std::endl << "* SURFACE / SURFACE" << std::endl;
  auto pair = bench_splitter<Kernel, Traits, Splitter>(
    verbose, num_iters, bucket_size, k, point_set, point_set, true);
  out.push_back(pair.first);
  out.push_back(pair.second);

  if (verbose) std::cout << std::endl << "* SURFACE / RANDOM" << std::endl;
  pair = bench_splitter<Kernel, Traits, Splitter>(
    verbose, num_iters, bucket_size, k, point_set, query_points);
  out.push_back(pair.first);
  out.push_back(pair.second);

  if (verbose) std::cout << std::endl << "* RANDOM / RANDOM" << std::endl;
  pair = bench_splitter<Kernel, Traits, Splitter>(
    verbose, num_iters, bucket_size, k, query_points, query_points);
  out.push_back(pair.first);
  out.push_back(pair.second);

  if (verbose) std::cout << std::endl;
  return out;
}

template<typename Kernel, typename Traits, typename Splitter>
std::vector<double> bench_translated(
  const std::string filename, const bool verbose,
  typename Kernel::FT d, std::size_t num_iters,
  std::size_t bucket_size, std::size_t k) {

  if (verbose) {
    std::cout << std::endl;
    std::cout << "---- TRANSLATION BENCH" << std::endl;
    std::cout << "- filename " << filename << " ... " << std::endl << std::endl;
  }

  std::vector<double> out;
  out.push_back(k);
  out.push_back(bucket_size);

  using Point_3 = typename Kernel::Point_3;
  using Vector_3 = typename Kernel::Vector_3;

  // Read point set.
  Timer timer;
  timer.start();
  std::vector<Point_3> point_set;
  std::ifstream in(filename);
  CGAL::IO::set_ascii_mode(in);
  Point_3 p;
  while (in >> p) {
    point_set.push_back(p);
  }
  assert(point_set.size() > 0);
  timer.stop();
  out.push_back(timer.time());

  if (verbose) {
    std::cout << "- PARAMETERS: " << std::endl;
    std::cout << "-            num k: " << k << std::endl;
    std::cout << "-      bucket size: " << bucket_size << std::endl;
    std::cout << "- num input points: " << point_set.size() << std::endl;
    std::cout << "- num query points: " << point_set.size() << std::endl;
    std::cout << std::endl;

    std::cout << "- CREATE POINT TIME: " << out.back() << " sec." << std::endl;

    std::ofstream outfile("0-original.xyz");
    CGAL::IO::write_XYZ(outfile, point_set);
  }

  // Translate point set.
  timer.reset();
  timer.start();
  auto query_points = point_set;
  const Vector_3 translation(d, d, d);
  translate_point_set<Kernel>(verbose, translation, query_points);
  CGAL_assertion(query_points.size() == point_set.size());
  timer.stop();
  out.push_back(timer.time());

  if (verbose) {
    std::cout << "- CREATE QUERY TIME: " << out.back() << " sec." << std::endl;
  }

  // Benching the tree.
  if (verbose) std::cout << std::endl << "* SURFACE / SURFACE" << std::endl;
  auto pair = bench_splitter<Kernel, Traits, Splitter>(
    verbose, num_iters, bucket_size, k, point_set, point_set, true);
  out.push_back(pair.first);
  out.push_back(pair.second);

  if (verbose) std::cout << std::endl << "* SURFACE / TRANSLATED" << std::endl;
  pair = bench_splitter<Kernel, Traits, Splitter>(
    verbose, num_iters, bucket_size, k, point_set, query_points);
  out.push_back(pair.first);
  out.push_back(pair.second);

  if (verbose) std::cout << std::endl << "* TRANSLATED / TRANSLATED" << std::endl;
  pair = bench_splitter<Kernel, Traits, Splitter>(
    verbose, num_iters, bucket_size, k, query_points, query_points);
  out.push_back(pair.first);
  out.push_back(pair.second);

  if (verbose) std::cout << std::endl;
  return out;
}

template<typename Kernel, typename Traits, typename Splitter>
void user_manual_bench(
  const std::string filename, const bool verbose, const std::size_t num_iters) {

  // Parameters:
  // std::size_t num_query_points; // number of query points; if 0 then it equals to the number of input points
  // std::size_t bucket_size;      // max bucket size per tree leaf
  // std::size_t k;                // number of k neighbors

  std::cout << std::endl << " --- USER MANUAL BENCH --- " << std::endl << std::endl;

  std::cout << "- Kernel: " << boost::typeindex::type_id<Kernel>() << std::endl;
  std::cout << "- Splitter: " << boost::typeindex::type_id<Splitter>() << std::endl;
  std::cout << std::endl;

  std::vector< std::vector<double> > out;
  if (!verbose) {
    std::cout << "k | ";
    std::cout << "bucket size | ";
    std::cout << "create points | ";
    std::cout << "create queris | ";
    std::cout << "surf/surf avg build | ";
    std::cout << "surf/surf avg search | ";
    std::cout << "surf/rndm avg build | ";
    std::cout << "surf/rndm avg search | ";
    std::cout << "rndm/rndm avg build | ";
    std::cout << "rndm/rndm avg search";
    std::cout << std::endl;
  }

  out.push_back(bench_random_in_bbox<Kernel, Traits, Splitter>(filename, verbose, 0, num_iters, 10, 10));
  out.push_back(bench_random_in_bbox<Kernel, Traits, Splitter>(filename, verbose, 0, num_iters, 20, 10));
  out.push_back(bench_random_in_bbox<Kernel, Traits, Splitter>(filename, verbose, 0, num_iters, 10, 20));
  out.push_back(bench_random_in_bbox<Kernel, Traits, Splitter>(filename, verbose, 0, num_iters, 20, 20));
  out.push_back(bench_random_in_bbox<Kernel, Traits, Splitter>(filename, verbose, 0, num_iters, 10, 30));
  out.push_back(bench_random_in_bbox<Kernel, Traits, Splitter>(filename, verbose, 0, num_iters, 20, 30));

  if (!verbose) {
    for (std::size_t k = 0; k < out.size(); ++k) {
      for (std::size_t i = 0; i < out[k].size(); ++i) std::cout << out[k][i] << " ";
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

template<typename Kernel, typename Traits, typename Splitter>
void translated_bench(
  const std::string filename, const bool verbose, const std::size_t num_iters) {

  // Parameters:
  // double d;                // offset distance between original point set and query points
  // std::size_t bucket_size; // max bucket size per tree leaf
  // std::size_t k;           // number of k neighbors

  std::cout << std::endl << " --- TRANSLATED BENCH --- " << std::endl << std::endl;

  std::cout << "- Kernel: " << boost::typeindex::type_id<Kernel>() << std::endl;
  std::cout << "- Splitter: " << boost::typeindex::type_id<Splitter>() << std::endl;
  std::cout << std::endl;

  std::vector< std::vector<double> > out;
  if (!verbose) {
    std::cout << "k | ";
    std::cout << "bucket size | ";
    std::cout << "create points | ";
    std::cout << "create queris | ";
    std::cout << "surf/surf avg build | ";
    std::cout << "surf/surf avg search | ";
    std::cout << "surf/trnl avg build | ";
    std::cout << "surf/trnl avg search | ";
    std::cout << "trnl/trnl avg build | ";
    std::cout << "trnl/trnl avg search";
    std::cout << std::endl;
  }

  out.push_back(bench_translated<Kernel, Traits, Splitter>(filename, verbose, 0.05, num_iters, 10, 10));
  out.push_back(bench_translated<Kernel, Traits, Splitter>(filename, verbose, 0.05, num_iters, 20, 10));
  out.push_back(bench_translated<Kernel, Traits, Splitter>(filename, verbose, 0.05, num_iters, 10, 20));
  out.push_back(bench_translated<Kernel, Traits, Splitter>(filename, verbose, 0.05, num_iters, 20, 20));
  out.push_back(bench_translated<Kernel, Traits, Splitter>(filename, verbose, 0.05, num_iters, 10, 30));
  out.push_back(bench_translated<Kernel, Traits, Splitter>(filename, verbose, 0.05, num_iters, 20, 30));

  if (!verbose) {
    for (std::size_t k = 0; k < out.size(); ++k) {
      for (std::size_t i = 0; i < out[k].size(); ++i) std::cout << out[k][i] << " ";
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

int main(int argc, char* argv[]) {

  // Parameters:
  std::cout.precision(20);
  const bool verbose = false; // should we export data
  const std::string filename = (argc > 1 ? argv[1] : "data/sphere.xyz"); // input file
  const std::size_t num_iters = 1; // number of iterations to get average timing

  user_manual_bench<SCD, SCD_traits, CGAL::Balanced_splitter<SCD_traits> >(
    filename, verbose, num_iters);

  // user_manual_bench<SCD, SCD_traits, CGAL::Sliding_midpoint<SCD_traits> >(
  //   filename, verbose, num_iters);

  translated_bench<SCD, SCD_traits, CGAL::Balanced_splitter<SCD_traits> >(
    filename, verbose, num_iters);

  // translated_bench<SCD, SCD_traits, CGAL::Sliding_midpoint<SCD_traits> >(
  //   filename, verbose, num_iters);

  return EXIT_SUCCESS;
}
