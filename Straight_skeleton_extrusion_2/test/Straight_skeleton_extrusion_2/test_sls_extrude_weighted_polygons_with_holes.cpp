#include "sls_extrusion_test_helpers.h"

#include <CGAL/extrude_skeleton.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/create_weighted_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Random.h>
#include <CGAL/random_polygon_2.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace SST = CGAL::Straight_skeletons_2::Extrusion_tests;

template <typename K>
using Case_list = std::vector<std::tuple<std::string, SST::Extrusion_input<K>, SST::Expectations> >;

// Weighted polygon with holes generated as in Straight_skeleton_2/test/Straight_skeleton_2/test_sls_weighted_polygons_with_holes.cpp:
// a long rectangle with fast vertical sides, and 2 * `hole_n` random holes with slow edges. The skeleton faces
// of the vertical sides surround the faces of the holes, so the builder adds artificial bisectors.
std::pair<CGAL::Polygon_with_holes_2<EPICK>, std::vector<std::vector<double> > >
generate_holed_rectangle(const int hole_n, const int hole_nv, const unsigned int seed)
{
  using Point_2 = EPICK::Point_2;
  using Polygon_2 = CGAL::Polygon_2<EPICK>;

  CGAL::Random rnd(seed);

  std::vector<Point_2> ob = { Point_2(-hole_n-1, -0.5), Point_2(hole_n+1, -0.5), Point_2(hole_n+1, 0.5), Point_2(-hole_n-1, 0.5) };
  CGAL::Polygon_with_holes_2<EPICK> pwh(Polygon_2(ob.begin(), ob.end()));

  std::vector<std::vector<double> > weights(1);
  weights[0].push_back(rnd.get_double(0.05, 0.5));
  weights[0].push_back(rnd.get_double(1, 10));
  if(rnd.get_int(0, 2) == 1)
    weights[0].push_back(rnd.get_double(0.05, 0.5));
  else
    weights[0].push_back(weights[0][0]);
  weights[0].push_back(rnd.get_double(1, 10));

  for(int hole_i=0; hole_i<hole_n; ++hole_i)
  {
    for(const int side : {1, -1})
    {
      Polygon_2 poly;
      CGAL::random_polygon_2(hole_nv, std::back_inserter(poly), CGAL::Random_points_in_square_2<Point_2>(0.25, rnd));
      std::vector<Point_2> hole_pts = poly.container();
      const double y_nudge = rnd.get_double(-0.25, 0.25);
      for(Point_2& pt : hole_pts)
        pt = Point_2(pt.x() + side * (hole_i + 1), pt.y() + y_nudge);
      // holes must be clockwise; same as Polygon_2::reverse_orientation(), which keeps the first vertex
      if(CGAL::polygon_area_2(hole_pts.begin(), hole_pts.end(), EPICK()) > 0)
        std::reverse(std::next(hole_pts.begin()), hole_pts.end());
      pwh.add_hole(Polygon_2(hole_pts.begin(), hole_pts.end()));
      weights.emplace_back(hole_pts.size(), 10.);
    }
  }

  return { pwh, weights };
}

template <typename K>
Case_list<K> weighted_polygons_with_holes_cases()
{
  using FT = typename K::FT;

  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  using Input = SST::Extrusion_input<K>;
  using Weights = typename Input::Weights;
  using Expectations = SST::Expectations;

  Case_list<K> cases;
  auto add = [&](const std::string& name, const Input& in, const Expectations& e)
  {
    cases.emplace_back(name, in, e);
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Non-simply connected skeleton faces

  // The fast right edge sweeps around the hole: its face surrounds the faces of the hole. The builder links them
  // with an artificial bisector, from (5.56, 1.06), t = 0.44 to (0.91, 1.06), t = 0.91, which the (single) boundary
  // of the face goes along in both directions. H = 0.3 is below the artificial bisector, H = 0.45 cuts it.
  for(double H : {0., 0.3, 0.45, -0.6})
  {
    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 4));
    pwh.add_hole(SST::make_rectangle<K>(6, 1.5, 7, 2.5, false /*cw*/));
    Input in(pwh);
    in.with_weights(Weights{{FT(1), FT(1), FT(0.1), FT(1)}, {FT(1), FT(1), FT(1), FT(1)}});
    if(H != 0)
      in.with_height(H);
    Expectations e;
    e.genus = 1;
    e.rel_tol = 1e-6;
    add("annular skeleton face, H = " + std::to_string(H), in, e);
  }

  // An artificial bisector that ends at an actual skeleton vertex, or very close to one. The fast right edge sweeps
  // around the hole (near the right edge), and the artificial bisector is shot from the leftmost vertex of the faces
  // of the hole along y = 2, until it hits the border of the face of the right edge:
  //  - w = 0.3: at the vertex (3.33, 2, t = 2) where the faces of the right, top, and bottom edges meet;
  //  - w = 0.25: at the vertex (2, 2, t = 2) where these faces and the face of the left edge meet;
  //  - w = 0.25 - 1e-12: on the bisector of the left and right edges, at about 6e-12 from (2, 2, t = 2).
  // H = 2 is the time of these vertices.
  for(const double w : {0.3, 0.25, 0.25 - 1e-12})
  {
    // Reference skeleton, to check that the configuration is as described
    {
      CGAL::Polygon_with_holes_2<EPICK> i_pwh(SST::make_rectangle<EPICK>(0, 0, 10, 4));
      i_pwh.add_hole(SST::make_polygon<EPICK>({{8.5, 2}, {9, 2.5}, {9.5, 2}, {9, 1.5}}));
      std::vector<std::vector<double> > i_weights = {{1, 1, w, 1}, {1, 1, 1, 1}};
      const auto ss = CGAL::create_interior_weighted_straight_skeleton_2(i_pwh, i_weights, EPICK());
      CGAL_SLS_CHECK_MSG(ss != nullptr, "w = " << w << ": no skeleton");
      if(!ss)
        continue;

      int artificial_bisectors = 0, to_actual_vertex = 0;
      for(auto h = ss->halfedges_begin(); h != ss->halfedges_end(); ++h)
      {
        if(h->face() != h->opposite()->face() || h->opposite() < h)
          continue;
        ++artificial_bisectors;
        // the end that is not the extreme vertex of the hole is the farthest from the contour edge
        auto v = (h->vertex()->time() > h->opposite()->vertex()->time()) ? h->vertex() : h->opposite()->vertex();
        if(v->trisegment()->e0().id() != v->trisegment()->e1().id()) // not the marker of artificial vertices
          ++to_actual_vertex;
      }
      CGAL_SLS_CHECK_MSG(artificial_bisectors == 1 && to_actual_vertex == ((w == 0.25 - 1e-12) ? 0 : 1),
                         "w = " << w << ": " << artificial_bisectors << " artificial bisector(s), "
                                << to_actual_vertex << " to an actual vertex");
    }

    for(const double H : {0., 0.5, 1., 2.})
    {
      Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 10, 4));
      pwh.add_hole(SST::make_polygon<K>({{8.5, 2}, {9, 2.5}, {9.5, 2}, {9, 1.5}})); // clockwise
      Input in(pwh);
      in.with_weights(Weights{{FT(1), FT(1), FT(w), FT(1)}, {FT(1), FT(1), FT(1), FT(1)}});
      if(H != 0)
        in.with_height(H);
      Expectations e;
      e.genus = 1;
      e.rel_tol = 1e-6;
      add("artificial bisector to an actual vertex (w = " + std::to_string(w) + "), H = " + std::to_string(H), in, e);
    }
  }

  // `generate_holed_rectangle()`, converted to K
  auto holed_rectangle = [](const int hole_n, const int hole_nv, const unsigned int seed)
  {
    const auto [i_pwh, i_weights] = generate_holed_rectangle(hole_n, hole_nv, seed);

    CGAL::Cartesian_converter<EPICK, K> to_k;
    Polygon_2 outer;
    for(const auto& p : i_pwh.outer_boundary())
      outer.push_back(to_k(p));
    Polygon_with_holes_2 pwh(outer);
    for(const auto& h : i_pwh.holes())
    {
      Polygon_2 hole;
      for(const auto& p : h)
        hole.push_back(to_k(p));
      pwh.add_hole(hole);
    }

    Weights weights;
    for(const auto& ws : i_weights)
      weights.emplace_back(ws.begin(), ws.end());

    return std::make_pair(pwh, weights);
  };

  // Several artificial bisectors in the same face: chained through holes (hole -> hole -> ... -> border of the face),
  // and several of them incident to the same hole. The face is then cut into an outer border and nested holes.
  //  - seed 33: nesting depth 3, and holes incident to three artificial bisectors
  //  - seed 47: nesting depth 6
  for(const unsigned int seed : {33u, 47u})
  {
    const auto [i_pwh, i_weights] = generate_holed_rectangle(4, 4, seed);

    // Reference skeleton, to check that the configuration is as described, and to get the maximal time
    std::vector<std::vector<EPICK::FT> > i_weights_ft(i_weights.begin(), i_weights.end());
    const auto ss = CGAL::create_interior_weighted_straight_skeleton_2(i_pwh, i_weights_ft, EPICK());
    CGAL_SLS_CHECK_MSG(ss != nullptr, "seed " << seed << ": no skeleton");
    if(!ss)
      continue;

    double t_max = 0;
    for(auto v = ss->vertices_begin(); v != ss->vertices_end(); ++v)
      t_max = (std::max)(t_max, CGAL::to_double(v->time()));

    const auto [pwh, weights] = holed_rectangle(4, 4, seed);

    // Heights are generic: they do not coincide with the time of a skeleton vertex
    for(const double f : {0., 0.25, 0.5, 0.75})
    {
      Input in(pwh);
      in.with_weights(weights);
      if(f != 0)
        in.with_height(f * t_max);
      Expectations e;
      e.genus = int(pwh.number_of_holes());
      e.rel_tol = 1e-6;
      add("nested artificial bisectors (seed " + std::to_string(seed) + "), " +
          ((f == 0) ? std::string("no height") : "H = " + std::to_string(f) + " * t_max"), in, e);
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Faces whose triangulation fails with EPICK

  // With an inexact kernel, rounded points can make the constrained triangulation of the border of a face fail:
  // the face is then triangulated again with the exact kernel.
  const std::vector<std::tuple<int, int, unsigned int, double> > epick_failures = {
    { 1, 4, 7u, 0.46571259654416619 },
    { 1, 4, 4u, 0.64419988942447493 },
    { 2, 4, 6u, 0.82955592923847854 },
    { 2, 4, 2u, 0.11994892904749498 } };
  for(const auto& [hole_n, hole_nv, seed, H] : epick_failures)
  {
    const auto [pwh, weights] = holed_rectangle(hole_n, hole_nv, seed);
    Input in(pwh);
    in.with_weights(weights).with_height(H);
    Expectations e;
    e.rel_tol = 1e-6;
    add("faces whose triangulation fails with EPICK (?) (" + std::to_string(hole_n) + ", " + std::to_string(hole_nv) +
        ", " + std::to_string(seed) + "), H = " + std::to_string(H), in, e);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Random holed rectangles

  // Heights: no height, fractions of the maximal time, and the times (computed with doubles) of the endpoints
  // of the artificial bisectors, where an artificial vertex or a hole vertex is (almost) on the offset.
  // With EPICK, the output can self-intersect at the latter heights.
  for(const auto& [hole_n, hole_nv] : std::vector<std::pair<int, int> >{ {1, 4}, {2, 4}, {1, 10} })
  {
    for(unsigned int seed=1; seed<=8; ++seed)
    {
      const auto [i_pwh, i_weights] = generate_holed_rectangle(hole_n, hole_nv, seed);
      std::vector<std::vector<EPICK::FT> > i_weights_ft(i_weights.begin(), i_weights.end());
      const auto ss = CGAL::create_interior_weighted_straight_skeleton_2(i_pwh, i_weights_ft, EPICK());
      if(!ss)
        continue;

      // Times of the endpoints of the artificial bisectors (in doubles), sorted and without duplicates
      std::set<double> times;
      for(auto h = ss->halfedges_begin(); h != ss->halfedges_end(); ++h)
        if(h->face() == h->opposite()->face())
          times.insert(CGAL::to_double(h->vertex()->time()));

      const std::vector<double> art_times(times.begin(), times.end());
      if(art_times.empty())
        continue;

      double t_max = 0;
      for(auto v = ss->vertices_begin(); v != ss->vertices_end(); ++v)
        t_max = (std::max)(t_max, CGAL::to_double(v->time()));

      std::vector<double> heights = { 0., 0.25 * t_max, 0.5 * t_max, 0.75 * t_max };
      heights.insert(heights.end(), art_times.begin(), art_times.end());

      const auto [pwh, weights] = holed_rectangle(hole_n, hole_nv, seed);
      for(const double H : heights)
      {
        Input in(pwh);
        in.with_weights(weights);
        if(H != 0)
          in.with_height(H);
        Expectations e;
        e.rel_tol = 1e-6;
        add("random holed rectangle (" + std::to_string(hole_n) + ", " + std::to_string(hole_nv) + ", " +
            std::to_string(seed) + "), " + ((H == 0) ? std::string("no height") : "H = " + std::to_string(H)), in, e);
      }
    }
  }

  return cases;
}

template <typename K>
std::vector<CGAL::Surface_mesh<typename K::Point_3> > run_cases()
{
  std::vector<CGAL::Surface_mesh<typename K::Point_3> > meshes;
  for(const auto& [name, in, e] : weighted_polygons_with_holes_cases<K>())
  {
    std::cout << "run " << name << " (" << typeid(K).name() << ")" << std::endl;
    meshes.emplace_back();
    const SST::Outcome o = extrude(in, meshes.back());
    CGAL_SLS_CHECK_MSG(o == SST::Outcome::SUCCESS, name + ": " + to_string(o));
    if(o == SST::Outcome::SUCCESS)
      check_extrusion(name, in, meshes.back(), e);
  }
  return meshes;
}

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  std::cout << "== EPICK ==" << std::endl;
  const auto epick_meshes = run_cases<EPICK>();
  std::cout << "== EPECK ==" << std::endl;
  const auto epeck_meshes = run_cases<EPECK>();

  // Cross-kernel consistency: same volume up to 1e-9 * input size
  std::cout << "== EPICK vs EPECK ==" << std::endl;
  const auto names = weighted_polygons_with_holes_cases<EPICK>();
  for(std::size_t i=0; i<epick_meshes.size(); ++i)
  {
    const std::string& name = std::get<0>(names[i]);
    if(num_faces(epick_meshes[i]) == 0 || num_faces(epeck_meshes[i]) == 0)
      continue; // already reported
    const double tol = 1e-9 * (std::max)(1e-300, SST::diag_xy(PMP::bbox(epeck_meshes[i])));
    CGAL_SLS_CHECK_MSG(SST::approx_equal(PMP::volume(epick_meshes[i]), PMP::volume(epeck_meshes[i]), tol),
                       "EPICK vs EPECK: " + name);
  }

  return SST::test_exit_code();
}
