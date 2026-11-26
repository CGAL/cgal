// Tests that the kd-tree query of hot-pixel snap rounding reports every hot
// pixel that a segment intersects. The exact test is applied to all hot pixels
// (brute force), and each intersected hot pixel must be among the candidates
// that the kd-tree returns.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <iterator>
#include <list>
#include <random>
#include <string>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Hot_pixel_snap_rounding_traits_2.h>
#include <CGAL/hot_pixel_snap_rounding_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits = CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>;
using FT = Kernel::FT;
using Point_2 = Kernel::Point_2;
using Segment_2 = Kernel::Segment_2;

namespace SR_internal = CGAL::Snap_Rounding_2::internal;

using Hot_pixel = SR_internal::Hot_pixel<Traits>;
using Segment_data = SR_internal::Segment_data<Traits>;
using Kd_tree = SR_internal::Hot_pixel_kd_tree<Traits, Hot_pixel*>;
using SEG_Direction = SR_internal::SEG_Direction;

// Queries the input segments and, as iterated snap rounding does, segments
// that connect hot-pixel centers.
bool test(const std::string& name, const std::vector<Segment_2>& segments, const FT& pixel_size,
          std::size_t num_center_segments, std::mt19937& gen) {
  Traits traits;
  std::list<Point_2> points;
  CGAL::compute_intersection_points(segments.begin(), segments.end(), std::back_inserter(points), true, traits);

  std::deque<Hot_pixel> pixels;
  std::vector<Hot_pixel*> hot_pixels;
  Kd_tree::Point_saved_pair_list pairs;
  for (const auto& p : points) {
    Hot_pixel* hp = &pixels.emplace_back(p, pixel_size, traits);
    hot_pixels.push_back(hp);
    pairs.emplace_back(hp->get_center(), hp);
  }
  Kd_tree tree(pairs, pixel_size);

  std::vector<Segment_2> queries(segments);
  std::uniform_int_distribution<std::size_t> pick(0, hot_pixels.size() - 1);
  for (std::size_t i = 0; i < num_center_segments; ++i) {
    Point_2 p = hot_pixels[pick(gen)]->get_center();
    Point_2 q = hot_pixels[pick(gen)]->get_center();
    if (p != q) queries.emplace_back(p, q);
  }

  std::size_t num_intersections = 0, num_misses = 0;
  std::vector<Hot_pixel*> candidates;
  for (const auto& s : queries) {
    Segment_data sd(s.source(), s.target());
    const SEG_Direction dir = sd.direction(traits);
    candidates.clear();
    tree.get_intersecting_points(candidates, s);
    std::sort(candidates.begin(), candidates.end());
    for (auto* hp : hot_pixels) {
      if (! hp->intersect(sd, dir, traits)) continue;
      ++num_intersections;
      if (std::binary_search(candidates.begin(), candidates.end(), hp)) continue;
      ++num_misses;
      std::cerr << name << ": the query of segment " << s << " missed the hot pixel centered at "
                << hp->get_center() << std::endl;
    }
  }

  std::cout << name << ": " << queries.size() << " queries, " << hot_pixels.size() << " hot pixels, "
            << num_intersections << " intersections, " << num_misses << " misses" << std::endl;
  return (num_misses == 0) && (num_intersections > 0);
}

// Segments on pixel boundaries and through pixel corners.
std::vector<Segment_2> boundary_segments() {
  std::vector<Segment_2> segs;
  const FT h(0.5);
  for (int k = -3; k <= 3; ++k) {
    segs.emplace_back(Point_2(k + h, -6), Point_2(k + h, 6));
    segs.emplace_back(Point_2(-6, k + h), Point_2(6, k + h));
    segs.emplace_back(Point_2(k - h, -h), Point_2(k + 5 + h, 5 + h));
    segs.emplace_back(Point_2(k + h, 5 + h), Point_2(k + 6 + h, -h));
  }
  return segs;
}

// Segments with endpoints on the half-integer grid; many degeneracies.
std::vector<Segment_2> grid_segments(std::size_t n, std::mt19937& gen) {
  std::uniform_int_distribution<int> d(-20, 20);
  std::vector<Segment_2> segs;
  for (std::size_t i = 0; i < n; ++i)
    segs.emplace_back(Point_2(FT(d(gen)) / 2, FT(d(gen)) / 2), Point_2(FT(d(gen)) / 2, FT(d(gen)) / 2));
  return segs;
}

// Long segments of a given direction and short cross segments, far from the origin.
std::vector<Segment_2> slanted_segments(double degrees, double offset, double scale, std::mt19937& gen) {
  const double a = degrees * CGAL_PI / 180.0;
  const double ux = std::cos(a), uy = std::sin(a);
  std::uniform_real_distribution<double> jitter(0, scale);
  std::vector<Segment_2> segs;
  for (int i = 0; i < 6; ++i) {
    double x = offset - i * uy * 3 * scale, y = offset + i * ux * 3 * scale;
    segs.emplace_back(Point_2(x + jitter(gen), y + jitter(gen)),
                      Point_2(x + 40 * scale * ux + jitter(gen), y + 40 * scale * uy + jitter(gen)));
  }
  for (int i = 0; i < 8; ++i) {
    double x = offset + 5 * i * scale * ux, y = offset + 5 * i * scale * uy;
    segs.emplace_back(Point_2(x + uy * scale, y - ux * scale), Point_2(x - 16 * uy * scale, y + 16 * ux * scale));
  }
  return segs;
}

int main() {
  std::mt19937 gen(17);
  bool ok = true;

  ok &= test("boundaries, pixel 1", boundary_segments(), FT(1), 30, gen);
  ok &= test("boundaries, pixel 1/3", boundary_segments(), FT(1) / FT(3), 30, gen);

  const FT pixel_sizes[] = {FT(1), FT(1) / FT(3), FT(0.1), FT(2.5)};
  for (const FT& pixel_size : pixel_sizes)
    ok &= test("half-integer grid, pixel " + std::to_string(CGAL::to_double(pixel_size)),
               grid_segments(30, gen), pixel_size, 30, gen);

  for (double degrees : {0.0, 35.0, 90.0, 119.0})
    ok &= test("slanted " + std::to_string(degrees), slanted_segments(degrees, 4e8, 3, gen), FT(1), 30, gen);

  // Rounding errors of the double approximations are comparable to the pixel size.
  ok &= test("large coordinates 1e9", slanted_segments(35, 1e9, 1e-5, gen), FT(1e-6), 30, gen);
  ok &= test("large coordinates 1e12", slanted_segments(119, 1e12, 1e-3, gen), FT(1e-3) / FT(3), 30, gen);

  if (! ok) {
    std::cerr << "FAILED" << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
