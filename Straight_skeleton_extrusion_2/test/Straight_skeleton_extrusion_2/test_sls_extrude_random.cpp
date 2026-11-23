#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Random.h>

#include "sls_extrusion_test_helpers.h"

#include <CGAL/extrude_skeleton.h>

#include <array>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace SSE = CGAL::Straight_skeleton_extrusion;
namespace SST = CGAL::Straight_skeletons_2::Extrusion_tests;

// Kernel-independent description of a random input, so that EPICK and EPECK see identical data
struct Random_input
{
  std::vector<std::vector<std::array<double, 2> > > contours; // [0] = outer (CCW), others = holes (CW)
  std::vector<std::vector<double> > weights;
  double height_fraction; // in ]0, 1.2[ of the unbounded height
  double outward_height;
};

Random_input random_input(CGAL::Random& rng, const bool grid)
{
  using Point_2 = EPICK::Point_2;
  Random_input ri;

  // outer boundary
  std::vector<Point_2> pts;
  const int n = rng.get_int(5, grid ? 14 : 24);
  if(grid)
  {
    // random_polygon_2 on integer points (unique points are required)
    std::vector<Point_2> cand;
    while(cand.size() < std::size_t(n))
    {
      Point_2 p(rng.get_int(0, 12), rng.get_int(0, 12));
      if(std::find(cand.begin(), cand.end(), p) == cand.end())
        cand.push_back(p);
    }
    CGAL::random_polygon_2(n, std::back_inserter(pts), cand.begin());
  }
  else
  {
    CGAL::Random_points_in_square_2<Point_2> gen(10., rng);
    CGAL::random_polygon_2(n, std::back_inserter(pts), gen);
  }

  CGAL::Polygon_2<EPICK> outer(pts.begin(), pts.end());
  if(outer.is_clockwise_oriented())
    outer.reverse_orientation();

  ri.contours.emplace_back();
  for(const Point_2& p : outer)
    ri.contours.back().push_back({p.x(), p.y()});

  // holes: random squares strictly inside the outer boundary, pairwise disjoint (with a margin)
  const CGAL::Bbox_2 bb = outer.bbox();
  const double size = (std::min)(bb.xmax() - bb.xmin(), bb.ymax() - bb.ymin());
  std::vector<EPICK::Iso_rectangle_2> holes;
  const int nh = rng.get_int(0, 4);
  for(int attempt=0; attempt<50 && int(holes.size())<nh; ++attempt)
  {
    double s = rng.get_double(0.05, 0.2) * size;
    double x = rng.get_double(bb.xmin(), bb.xmax() - s), y = rng.get_double(bb.ymin(), bb.ymax() - s);
    if(grid) {
      s = (std::max)(1., std::floor(s));
      x = std::floor(x);
      y = std::floor(y);
    }

    const EPICK::Iso_rectangle_2 r(Point_2(x, y), Point_2(x + s, y + s));
    const EPICK::Iso_rectangle_2 margin(Point_2(x - 0.1 * s, y - 0.1 * s), Point_2(x + 1.1 * s, y + 1.1 * s));

    bool ok = true;
    for(int i=0; i<4 && ok; ++i)
      ok = (outer.bounded_side(margin.vertex(i)) == CGAL::ON_BOUNDED_SIDE);
    for(auto e = outer.edges_begin(); e != outer.edges_end() && ok; ++e)
      ok = !CGAL::do_intersect(*e, margin);
    for(const auto& h : holes)
      ok = ok && !CGAL::do_intersect(h, margin);
    if(!ok)
      continue;

    holes.push_back(r);
    ri.contours.push_back({{x, y}, {x, y + s}, {x + s, y + s}, {x + s, y}}); // clockwise
  }

  for(const auto& c : ri.contours)
  {
    ri.weights.emplace_back();
    for(std::size_t i=0; i<c.size(); ++i)
      ri.weights.back().push_back(grid ? double(rng.get_int(1, 4)) : rng.get_double(1., 4.));

    // collinear consecutive edges cannot have the same weights
    const std::size_t m = c.size();
    for(int pass=0; pass<2; ++pass)
    {
      for(std::size_t k=0; k<m; ++k)
      {
        const auto& a = c[(k + m - 1) % m]; const auto& b = c[k]; const auto& d = c[(k + 1) % m];
        if(CGAL::collinear(Point_2(a[0], a[1]), Point_2(b[0], b[1]), Point_2(d[0], d[1])))
          ri.weights.back()[(k + 1) % m] = ri.weights.back()[k];
      }
    }
  }
  ri.height_fraction = rng.get_double(0.05, 1.2);
  ri.outward_height = rng.get_double(0.05, 3.);
  return ri;
}

template <typename K>
SST::Extrusion_input<K> to_input(const Random_input& ri,
                                 const double weight_sign)
{
  SST::Extrusion_input<K> in;
  auto contour = [](const std::vector<std::array<double, 2> >& c)
  {
    CGAL::Polygon_2<K> p;
    for(const auto& q : c)
      p.push_back(typename K::Point_2(q[0], q[1]));
    return p;
  };

  in.pwh = CGAL::Polygon_with_holes_2<K>(contour(ri.contours[0]));
  for(std::size_t i=1; i<ri.contours.size(); ++i)
    in.pwh.add_hole(contour(ri.contours[i]));

  for(std::size_t c=0; c<ri.weights.size(); ++c)
  {
    in.weights.emplace_back();
    for(double w : ri.weights[c])
      in.weights.back().push_back(typename K::FT(weight_sign * w));
  }

  return in;
}

template <typename K>
std::vector<CGAL::Surface_mesh<typename K::Point_3> > run_random_input(const std::string& name,
                                                                       const Random_input& ri)
{
  std::cout << "== run " << name << " ==" << std::endl;

  using Mesh = CGAL::Surface_mesh<typename K::Point_3>;
  std::vector<Mesh> meshes(4);
  const int genus = int(ri.contours.size()) - 1;

  // inward, unbounded
  SST::Extrusion_input<K> in = to_input<K>(ri, 1.);
  SST::Outcome o = extrude(in, meshes[0]);
  CGAL_SLS_CHECK_MSG(o == SST::Outcome::SUCCESS, name << " inward: " << to_string(o));
  if(o != SST::Outcome::SUCCESS)
    return meshes;
  SST::Expectations e;
  e.genus = genus;
  e.rel_tol = 1e-6;
  check_extrusion(name + " inward", in, meshes[0], e);
  const double full_height = PMP::bbox(meshes[0]).zmax();

  // inward, capped
  SST::Extrusion_input<K> inc = in;
  inc.with_height(ri.height_fraction * full_height);
  o = extrude(inc, meshes[1]);
  CGAL_SLS_CHECK_MSG(o == SST::Outcome::SUCCESS, name << " capped: " << to_string(o));
  if(o == SST::Outcome::SUCCESS)
  {
    SST::Expectations ec;
    ec.genus = genus;
    ec.rel_tol = 1e-6;
    check_extrusion(name + " capped", inc, meshes[1], ec);
    const double H = CGAL::to_double(*inc.height);
    CGAL_SLS_CHECK_MSG(SST::approx_equal(PMP::volume(meshes[1]), SST::clipped_volume<K>(meshes[0], H), 1e-6), name << " capped volume");
  }

  // downward: mirror of capped
  SST::Extrusion_input<K> ind = inc;
  ind.with_height(- *inc.height);
  o = extrude(ind, meshes[2]);
  CGAL_SLS_CHECK_MSG(o == SST::Outcome::SUCCESS, name << " downward: " << to_string(o));
  if(o == SST::Outcome::SUCCESS)
  {
    const double tol = 1e-9 * SST::diag_xy(PMP::bbox(meshes[1]));
    CGAL_SLS_CHECK_MSG(SST::approx_equal(PMP::volume(meshes[1]), PMP::volume(meshes[2]), tol),
                       name + " mirror");
  }

  // outward
  SST::Extrusion_input<K> out = to_input<K>(ri, -1.);
  out.with_height(ri.outward_height);
  o = extrude(out, meshes[3]);
  CGAL_SLS_CHECK_MSG(o == SST::Outcome::SUCCESS, name << " outward: " << to_string(o));
  if(o == SST::Outcome::SUCCESS)
  {
    SST::Expectations eo;
    eo.rel_tol = 1e-6;
    check_extrusion(name + " outward", out, meshes[3], eo);
  }
  return meshes;
}

// Usage: test_sls_extrude_random [seed] [number of polygons]
int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  const unsigned int seed = (argc > 1) ? static_cast<unsigned int>(std::atoi(argv[1])) : 20261112u;
  const int n = (argc > 2) ? std::atoi(argv[2]) : 60;
  std::cout << "Seed: " << seed << ", " << n << " polygons per family" << std::endl;
  CGAL::Random rng(seed);

  // generic: coordinates are arbitrary doubles;
  // grid: coordinates are integers in a small range (many parallel / collinear edges, simultaneous events).
  for(const bool grid : {false, true})
  {
    for(int i=0; i<n; ++i)
    {
      const Random_input ri = random_input(rng, grid);
      const std::string name = std::string(grid ? "grid" : "generic") + " #" + std::to_string(i) +
                               " (seed " + std::to_string(seed) + ")";
      const int before = SST::failure_count();

      const auto epick_meshes = run_random_input<EPICK>(name + " EPICK", ri);

      // every 4th input also with EPECK (slower), compared with EPICK
      if(i % 4 == 0)
      {
        const auto epeck_meshes = run_random_input<EPECK>(name + " EPECK", ri);
        for(std::size_t k=0; k<epick_meshes.size(); ++k)
        {
          if(CGAL::is_empty(epick_meshes[k]) || CGAL::is_empty(epeck_meshes[k])) // failure already reported
            continue;

          const double tol = 1e-9 * SST::diag_xy(PMP::bbox(epeck_meshes[k]));
          CGAL_SLS_CHECK_MSG(SST::approx_equal(PMP::volume(epick_meshes[k]), PMP::volume(epeck_meshes[k]), tol),
                             name + " EPICK vs EPECK, config " + std::to_string(k));
        }
      }

      if(SST::failure_count() != before) // print the input to reproduce the failure
      {
        std::cerr << "Failing input (" << name << "):" << std::endl;
        for(std::size_t c=0; c<ri.contours.size(); ++c)
        {
          std::cerr << (c == 0 ? "  outer:" : "  hole:");
          for(std::size_t j=0; j<ri.contours[c].size(); ++j)
            std::cerr << " (" << ri.contours[c][j][0] << ", " << ri.contours[c][j][1] << ") w=" << ri.weights[c][j];
          std::cerr << std::endl;
        }
        std::cerr << "  height fraction " << ri.height_fraction << ", outward height " << ri.outward_height << std::endl;
      }
    }
  }

  return SST::test_exit_code();
}
