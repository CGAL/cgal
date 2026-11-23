#include "sls_extrusion_test_helpers.h"

#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/extrude_skeleton.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <array>
#include <iostream>
#include <string>
#include <vector>

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace SSE = CGAL::Straight_skeleton_extrusion;
namespace SST = CGAL::Straight_skeletons_2::Extrusion_tests;

template <typename K>
struct Snapping_tests
{
  using FT = typename K::FT;
  using Point_2 = typename K::Point_2;
  using Point_3 = typename K::Point_3;

  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

  using Mesh = CGAL::Surface_mesh<Point_3>;
  using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;

  using Input = SST::Extrusion_input<K>;
  using Weights = typename Input::Weights;
  using Expectations = SST::Expectations;
  using Outcome = SST::Outcome;

  static constexpr bool is_exact = SSE::internal::is_exact_kernel<K>;

  static Mesh run(const std::string& name,
                  const Input& in,
                  const Expectations& e = Expectations())
  {
    std::cout << "run " << name << std::endl;

    Mesh m;
    const Outcome o = extrude(in, m);
    CGAL_SLS_CHECK_MSG(o == Outcome::SUCCESS, name + ": " + to_string(o));
    if(o == Outcome::SUCCESS)
      check_extrusion(name, in, m, e);
    return m;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Vertical slabs

  static double scale(const Mesh& m) { return SST::diag_xy(PMP::bbox(m)); }

  // signed distance (times |b-a|) of the xy-projection of `p` to the line (a, b)
  static FT line_side(const Point_2& a, const Point_2& b, const Point_3& p)
  {
    return (b.x() - a.x()) * (p.y() - a.y()) - (b.y() - a.y()) * (p.x() - a.x());
  }

  // Every vertex of `m` that is close to the vertical plane above the segment [a;b]
  // (i.e., whose projection is close to the segment) must be exactly on that plane:
  // these are the points of the vertical face, which have been snapped.
  //
  // For axis-aligned walls, the coordinate must be exactly the input coordinate, whatever the kernel.
  // Otherwise, it is exact for exact kernels, and up to rounding for inexact kernels.
  //
  // Returns the number of vertices on the wall.
  static std::size_t check_wall(const std::string& name,
                                const Mesh& m,
                                const Point_2& a, const Point_2& b)
  {
    const double len = CGAL::to_double(CGAL::approximate_sqrt(CGAL::squared_distance(a, b)));
    const double near_tol = 1e-6 * scale(m); // the slope of the approximated vertical faces is ~1e-9
    const double round_tol = 1e-13 * scale(m);

    std::size_t n = 0;
    for(vertex_descriptor v : vertices(m))
    {
      const Point_3& p = m.point(v);

      // the projection of p onto the line must be within the segment
      const double s = CGAL::to_double((p.x() - a.x()) * (b.x() - a.x()) + (p.y() - a.y()) * (b.y() - a.y())) / (len * len);
      if(s < -1e-6 || s > 1 + 1e-6)
        continue;

      const double d = std::abs(CGAL::to_double(line_side(a, b, p))) / len;
      if(d > near_tol)
        continue;

      ++n;
      if(a.x() == b.x())
        CGAL_SLS_CHECK_MSG(p.x() == a.x(), name << ": " << p << " should be on the vertical plane x = " << a.x());
      else if(a.y() == b.y())
        CGAL_SLS_CHECK_MSG(p.y() == a.y(), name << ": " << p << " should be on the vertical plane y = " << a.y());
      else if(is_exact)
        CGAL_SLS_CHECK_MSG(CGAL::is_zero(line_side(a, b, p)), name << ": " << p << " should be exactly above (" << a << ") (" << b << ")");
      else
        CGAL_SLS_CHECK_MSG(d <= round_tol, name << ": " << p << " is at " << d << " of the plane above (" << a << ") (" << b << ")");
    }
    return n;
  }

  static Point_2 P(double x, double y) { return Point_2(x, y); }

  // A single axis-aligned vertical wall: square [0,2]^2, the left edge (v3, v0) has weight 0.
  // The other fronts move at speed 1, so the cross-section at height z is [0, 2-z] x [z, 2-z]:
  // the ridge goes from (0,1,1) (on the wall) to (1,1,1).
  static void test_axis_aligned_wall()
  {
    const Polygon_2 sq = SST::make_rectangle<K>(0, 0, 2, 2);
    const Weights w{{FT(0), FT(1), FT(1), FT(1)}};
    auto volume = [](double H) { return 4 * H - 3 * H * H + 2 * H * H * H / 3; }; // integral of (2-z)(2-2z)

    {
      Input in(sq);
      in.with_weights(w);
      Expectations e;
      e.nv = 6;
      e.nf = 8;
      e.volume = volume(1);
      e.rel_tol = 1e-6;
      const Mesh m = run("axis-aligned wall, unbounded", in, e);
      CGAL_SLS_CHECK_MSG(check_wall("axis-aligned wall, unbounded", m, P(0, 2), P(0, 0)) == 3, "wall: ridge end + 2 corners");
      CGAL_SLS_CHECK_MSG(SST::has_vertex_approx(m, Point_3(0, 1, 1), 1e-6), "axis-aligned wall, unbounded: ridge end");
    }

    for(double H : {0.5, -0.5, 0.1, 1. - 1e-9})
    {
      const std::string name = "axis-aligned wall, H = " + std::to_string(H);
      Input in(sq);
      in.with_weights(w).with_height(H);
      Expectations e;
      e.nv = 8;
      e.nf = 12;
      e.volume = volume(std::abs(H));
      e.top_area = (2 - std::abs(H)) * (2 - 2 * std::abs(H));
      e.zmax_abs = std::abs(H);
      e.rel_tol = 1e-6;
      const Mesh m = run(name, in, e);
      CGAL_SLS_CHECK_MSG(check_wall(name, m, P(0, 2), P(0, 0)) == 4, name << ": wall: 2 top + 2 corners");

      // the top points on the wall: x is exact, z is exact
      for(vertex_descriptor v : vertices(m))
      {
        const Point_3& p = m.point(v);
        if(p.z() != 0)
          CGAL_SLS_CHECK_MSG(p.z() == FT(H), name << ": " << p << " should be at the height");
      }
    }
  }

  // A vertical wall whose edge is not axis-aligned: square of side 5 with 3-4-5 edges.
  // The opposite edge moves towards the wall at speed 1, the two others shrink the length along the wall:
  // the cross-section at height z has area (5 - z)(5 - 2z).
  static void test_rotated_wall()
  {
    const Polygon_2 sq = SST::make_polygon<K>({{0, 0}, {4, 3}, {1, 7}, {-3, 4}});
    const Weights w{{FT(0), FT(1), FT(1), FT(1)}}; // the edge (-3,4) -> (0,0)

    for(double H : {1., 2.25, -1.})
    {
      const std::string name = "rotated wall, H = " + std::to_string(H);
      const double h = std::abs(H);
      Input in(sq);
      in.with_weights(w).with_height(H);
      Expectations e;
      e.nv = 8;
      e.nf = 12;
      e.volume = 25 * h - 7.5 * h * h + 2 * h * h * h / 3;
      e.top_area = (5 - h) * (5 - 2 * h);
      e.rel_tol = 1e-6;
      const Mesh m = run(name, in, e);
      CGAL_SLS_CHECK_MSG(check_wall(name, m, P(-3, 4), P(0, 0)) == 4, name << ": wall: 2 top + 2 corners");
    }

    {
      Input in(sq);
      in.with_weights(w);
      Expectations e;
      e.nv = 6;
      e.volume = 25 * 2.5 - 7.5 * 2.5 * 2.5 + 2 * 2.5 * 2.5 * 2.5 / 3;
      e.rel_tol = 1e-6;
      const Mesh m = run("rotated wall, unbounded", in, e);
      CGAL_SLS_CHECK_MSG(check_wall("rotated wall, unbounded", m, P(-3, 4), P(0, 0)) == 3, "rotated wall, unbounded: ridge end + 2 corners");
    }
  }

  // Two adjacent vertical walls: square [0,2]^2, edges x = 0 and y = 0 are vertical.
  // Points on both walls are snapped onto their common contour vertex (0, 0).
  static void test_adjacent_walls()
  {
    const Polygon_2 sq = SST::make_rectangle<K>(0, 0, 2, 2);
    const Weights w{{FT(0), FT(0), FT(1), FT(1)}};

    for(double H : {1., -1., 0.3})
    {
      const std::string name = "adjacent walls, H = " + std::to_string(H);
      const double h = std::abs(H);
      Input in(sq);
      in.with_weights(w).with_height(H);
      Expectations e;
      e.nv = 8;
      e.nf = 12;
      e.volume = (8 - (2 - h) * (2 - h) * (2 - h)) / 3; // integral of (2-z)^2
      e.top_area = (2 - h) * (2 - h);
      e.rel_tol = 1e-6;
      const Mesh m = run(name, in, e);
      CGAL_SLS_CHECK_MSG(check_wall(name, m, P(0, 2), P(0, 0)) == 4, name << ": left wall");
      CGAL_SLS_CHECK_MSG(check_wall(name, m, P(0, 0), P(2, 0)) == 4, name << ": bottom wall");
      CGAL_SLS_CHECK_MSG(SST::has_vertex(m, Point_3(0, 0, H)), name << ": corner of both walls at the top");
    }

    // unbounded: a single apex, on both walls (the finite weight of the walls makes it slightly lower than 2)
    {
      Input in(sq);
      in.with_weights(w);
      Expectations e;
      e.nv = 5;
      e.nf = 6;
      e.volume = 8. / 3.;
      e.zmax_abs = 2;
      e.rel_tol = 1e-6;
      const Mesh m = run("adjacent walls, unbounded", in, e);
      bool found = false;
      for(vertex_descriptor v : vertices(m))
      {
        const Point_3& p = m.point(v);
        if(p.z() != 0)
        {
          found = true;
          CGAL_SLS_CHECK_MSG(p.x() == 0 && p.y() == 0, "adjacent walls, unbounded: apex " << p << " should be above (0,0)");
        }
      }
      CGAL_SLS_CHECK_MSG(found, "adjacent walls, unbounded: no apex");
    }
  }

  // Two non-adjacent vertical walls: L-shape with a chamfered reflex corner
  // (4,2) -> (2.5,2) [wall] -> (2,2.5) [chamfer] -> (2,4) [wall].
  // The face of the chamfer is a triangle whose apex is above the intersection (2,2) of the supporting lines
  // of the walls: it is snapped with the line intersection rather than with a common contour vertex.
  static void test_non_adjacent_walls()
  {
    const Polygon_2 L = SST::make_polygon<K>({{0, 0}, {4, 0}, {4, 2}, {2.5, 2}, {2, 2.5}, {2, 4}, {0, 4}});
    const Weights w{{FT(1), FT(1), FT(1), FT(0), FT(1), FT(0), FT(1)}};

    for(double H : {0.25, 1., -1., 0.})
    {
      const std::string name = "non-adjacent walls, H = " + std::to_string(H);
      Input in(L);
      in.with_weights(w);
      if(H != 0)
        in.with_height(H);
      Expectations e;
      e.rel_tol = 1e-6;
      const Mesh m = run(name, in, e);
      check_wall(name, m, P(4, 2), P(2.5, 2));
      check_wall(name, m, P(2, 2.5), P(2, 4));

      // the apex of the chamfer face, at height sqrt(2)/4, is above (2,2)
      if(H == 0 || std::abs(H) > 0.36)
      {
        bool found = false;
        for(vertex_descriptor v : vertices(m))
        {
          const Point_3& p = m.point(v);
          if(CGAL::abs(p.z()) > 0.35 && CGAL::abs(p.z()) < 0.36)
          {
            found = true;
            CGAL_SLS_CHECK_MSG(p.x() == 2 && p.y() == 2, name << ": chamfer apex " << p << " should be above (2,2)");
          }
        }
        CGAL_SLS_CHECK_MSG(found, name << ": chamfer apex not found");
      }
    }
  }

  // Vertical walls on a hole: a vertical shaft through a pyramid frustum.
  static void test_hole_walls()
  {
    Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 6, 6));
    pwh.add_hole(SST::make_rectangle<K>(2, 2, 4, 4, false /*cw*/));
    const Weights w{{FT(1), FT(1), FT(1), FT(1)}, {FT(0), FT(0), FT(0), FT(0)}};

    for(double H : {1., -1.5})
    {
      const std::string name = "vertical hole, H = " + std::to_string(H);
      const double h = std::abs(H);
      Input in(pwh);
      in.with_weights(w).with_height(H);
      Expectations e;
      e.nv = 16;
      e.genus = 1;
      e.volume = (216 - (6 - 2 * h) * (6 - 2 * h) * (6 - 2 * h)) / 6 - 4 * h; // integral of (6-2z)^2 - 4
      e.top_area = (6 - 2 * h) * (6 - 2 * h) - 4;
      e.rel_tol = 1e-6;
      const Mesh m = run(name, in, e);
      for(const auto& c : std::vector<std::array<double, 2> >{{2, 2}, {2, 4}, {4, 4}, {4, 2}})
        CGAL_SLS_CHECK_MSG(SST::has_vertex(m, Point_3(c[0], c[1], H)), name << ": top corner of the shaft above (" << c[0] << " " << c[1] << ")");
      CGAL_SLS_CHECK_MSG(check_wall(name, m, P(2, 2), P(2, 4)) == 4, name);
      CGAL_SLS_CHECK_MSG(check_wall(name, m, P(4, 4), P(4, 2)) == 4, name);
    }
  }

  // Outward extrusion with a vertical wall: the other fronts move away from the polygon.
  static void test_outward_wall()
  {
    const Polygon_2 sq = SST::make_rectangle<K>(0, 0, 2, 2);
    const Weights w{{FT(0), FT(-1), FT(-1), FT(-1)}};

    for(double H : {0.5, -0.5, 2.})
    {
      const std::string name = "outward wall, H = " + std::to_string(H);
      const double h = std::abs(H);
      Input in(sq);
      in.with_weights(w).with_height(H);
      Expectations e;
      e.nv = 8;
      e.nf = 12;
      e.volume = 4 * h + 3 * h * h + 2 * h * h * h / 3; // integral of (2+z)(2+2z)
      e.top_area = (2 + h) * (2 + 2 * h);
      e.rel_tol = 1e-6;
      const Mesh m = run(name, in, e);
      CGAL_SLS_CHECK_MSG(check_wall(name, m, P(0, 2 + h), P(0, -h)) == 4, name << ": wall: 2 top + 2 corners");
    }
  }

  // All edges vertical: the result must be exactly the prism, whatever the kernel.
  static void test_prisms()
  {
    auto check_prism = [](const std::string& name, const Mesh& m,
                          const Polygon_with_holes_2& pwh, const double H)
    {
      std::size_t n = pwh.outer_boundary().size();
      for(const auto& h : pwh.holes())
        n += h.size();
      CGAL_SLS_CHECK_MSG(num_vertices(m) == 2 * n, name << ": #V = " << num_vertices(m));

      for(vertex_descriptor v : vertices(m))
      {
        const Point_3& p = m.point(v);
        bool is_contour_vertex = false;
        auto test = [&](const Polygon_2& c) {
          for(const Point_2& q : c)
            if(q.x() == p.x() && q.y() == p.y())
              is_contour_vertex = true;
        };
        test(pwh.outer_boundary());
        for(const auto& h : pwh.holes())
          test(h);

        CGAL_SLS_CHECK_MSG(is_contour_vertex, name << ": " << p << " is not above a contour vertex");
        CGAL_SLS_CHECK_MSG(p.z() == 0 || p.z() == FT(H), name << ": " << p << " is neither at the bottom nor at the top");
      }
    };

    {
      const Polygon_2 r = SST::make_rectangle<K>(0, 0, 3, 1);
      for(double H : {1., -2.5, 0.1, 1e3})
      {
        const std::string name = "vertical rectangle, H = " + std::to_string(H);
        Input in(r);
        in.with_angles(Weights{{FT(90), FT(90), FT(90), FT(90)}}).with_height(H);
        Expectations e;
        e.nv = 8;
        e.nf = 12;
        e.volume = 3 * std::abs(H);
        e.top_area = 3;
        const Mesh m = run(name, in, e);
        check_prism(name, m, Polygon_with_holes_2(r), H);
      }
    }

    {
      Polygon_with_holes_2 pwh(SST::make_polygon<K>({{0, 0}, {5, 0}, {5, 2}, {2, 2}, {2, 5}, {0, 5}}));
      pwh.add_hole(SST::make_rectangle<K>(0.5, 0.5, 1.5, 1.5, false /*cw*/));
      const double area = 5 * 2 + 2 * 3 - 1;
      for(double H : {0.75, -3.})
      {
        const std::string name = "vertical L with hole, H = " + std::to_string(H);
        Input in(pwh);
        in.with_weights(Input::uniform(pwh, 0.)).with_height(H);
        Expectations e;
        e.genus = 1;
        e.volume = area * std::abs(H);
        e.top_area = area;
        const Mesh m = run(name, in, e);
        check_prism(name, m, pwh, H);
      }
    }

    // A rotated, non-convex polygon: the corners are snapped onto the contour vertices, and the projections onto
    // non-axis aligned walls must not move these points.
    {
      const Polygon_2 p = SST::make_polygon<K>({{0, 0}, {4, 3}, {2.5, 5}, {1, 4}, {-1.5, 4.5}, {-3, 4}});
      for(double H : {2., -0.5})
      {
        const std::string name = "vertical rotated polygon, H = " + std::to_string(H);
        Input in(p);
        in.with_weights(Input::uniform(Polygon_with_holes_2(p), 0.)).with_height(H);
        Expectations e;
        e.volume = CGAL::to_double(p.area()) * std::abs(H);
        e.top_area = CGAL::to_double(p.area());
        const Mesh m = run(name, in, e);
        check_prism(name, m, Polygon_with_holes_2(p), H);
      }
    }
  }

#if 0
  // A height so large that the faces of the vertical walls (which have a very large, but finite, weight) meet,
  // and points are on vertical slabs that have no common vertical line:
  // - all edges vertical: the ridge ends are on three walls (two of which are parallel), and are snapped
  //   with the first two (a warning is emitted);
  // - steep (but not vertical) short edges: the long vertical walls meet before the fronts of the short edges
  //   arrive, and the ridge ends are on two parallel walls: they are not snapped (a warning is emitted).
  static void test_height_beyond_vertical_slabs()
  {
    {
      const std::string name = "vertical rectangle, H = 1e12";
      Input in(SST::make_rectangle<K>(0, 0, 2, 1));
      in.with_angles(Weights{{FT(90), FT(90), FT(90), FT(90)}}).with_height(1e12);
      Expectations e;
      e.nv = 6;
      e.check_self_intersections = false; // very large coordinates
      run(name, in, e);
    }

    {
      // the angle i is for the edge (v_{i-1}, v_i): 0 and 2 are the short edges
      const std::string name = "parallel vertical walls, steep short edges";
      Input in(SST::make_rectangle<K>(0, 0, 1, 1e-4));
      in.with_angles(Weights{{FT(89.9999), FT(90), FT(89.9999), FT(90)}});
      Expectations e;
      e.nv = 6;
      run(name, in, e);
    }
  }
#endif

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Height clamping

  // The nodes of the (unweighted) skeleton of `pwh`, computed with an exact kernel
  static std::vector<std::pair<EPECK::Point_2, EPECK::FT> > skeleton_nodes(const Polygon_with_holes_2& pwh)
  {
    CGAL::Cartesian_converter<K, EPECK> to_exact;
    CGAL::Polygon_2<EPECK> outer;
    for(const Point_2& p : pwh.outer_boundary())
      outer.push_back(to_exact(p));
    std::vector<CGAL::Polygon_2<EPECK> > holes;
    for(const auto& h : pwh.holes())
    {
      holes.emplace_back();
      for(const Point_2& p : h)
        holes.back().push_back(to_exact(p));
    }

    auto ss = CGAL::create_interior_straight_skeleton_2(outer.vertices_begin(), outer.vertices_end(),
                                                        holes.begin(), holes.end(), EPECK());
    std::vector<std::pair<EPECK::Point_2, EPECK::FT> > nodes;
    for(auto v = ss->vertices_begin(); v != ss->vertices_end(); ++v)
      if(!v->is_contour())
        nodes.emplace_back(v->point(), v->time());
    return nodes;
  }

  // The nodes whose exact time is exactly `H`
  static std::vector<EPECK::Point_2> nodes_at_time(const Polygon_with_holes_2& pwh, const double H)
  {
    std::vector<EPECK::Point_2> nodes;
    for(const auto& n : skeleton_nodes(pwh))
      if(n.second == EPECK::FT(H))
        nodes.push_back(n.first);
    return nodes;
  }

  // Extrudes `pwh` at the height of skeleton nodes, for which the time computed with doubles might be
  // slightly off. Nothing can be above the height. If the nodes are on the top contour, they are exactly at the
  // height; otherwise (degenerate top contour), they can be slightly below.
  static void check_nodes_at_height(const std::string& name,
                                    const Polygon_with_holes_2& pwh,
                                    const double H,
                                    const std::size_t expected_nodes_n,
                                    const bool on_top_contour)
  {
    const std::vector<EPECK::Point_2> nodes = nodes_at_time(pwh, H);
    CGAL_SLS_CHECK_MSG(nodes.size() == expected_nodes_n,
                       name << ": " << nodes.size() << " nodes at the height, expected " << expected_nodes_n);

    for(const double sign : {1., -1.})
    {
      const std::string sname = name + ", H = " + std::to_string(sign * H);
      Input in(pwh);
      in.with_height(sign * H);
      Expectations e;
      e.zmax_abs = H;
      e.rel_tol = 1e-6;
      const Mesh m = run(sname, in, e);

      for(const EPECK::Point_2& n : nodes)
      {
        bool found = false;
        for(vertex_descriptor v : vertices(m))
        {
          const Point_3& p = m.point(v);
          if(SST::approx_equal(p.x(), n.x(), 1e-9, 1e-9) && SST::approx_equal(p.y(), n.y(), 1e-9, 1e-9))
          {
            if(p.z() == 0)
              continue;
            found = true;
            if(on_top_contour || is_exact)
              CGAL_SLS_CHECK_MSG(p.z() == FT(sign * H), sname << ": node " << p << " should be exactly at the height");
            else
              CGAL_SLS_CHECK_MSG(SST::approx_equal(p.z(), sign * H, 1e-12, 0),
                                 sname << ": node " << p << " should be at the height, up to rounding");
          }
        }
        CGAL_SLS_CHECK_MSG(found, sname << ": node (" << n << ") not found");
      }

      // nothing just below the height on the top contour: all points are either exactly at the height,
      // or at a distance that is not a numerical error
      if(on_top_contour)
      {
        for(vertex_descriptor v : vertices(m))
        {
          const double z = std::abs(CGAL::to_double(m.point(v).z()));
          CGAL_SLS_CHECK_MSG(z == H || z < H * (1 - 1e-9), sname << ": " << m.point(v) << " is just below the height");
        }
      }
    }
  }

  static void test_nodes_at_height()
  {
    // Apexes (degenerate top contour):
    // - 5-12-13 triangle: the time of the incenter computed with doubles is 1.9999999999999998, its exact time is 2
    // - 3-4-5 triangle: the time of the incenter computed with doubles is 1.0000000000000002, its exact time is 1
    // - 8-15-17 triangle: 3.0000000000000004 vs 3
    check_nodes_at_height("5-12-13 triangle", Polygon_with_holes_2(SST::make_polygon<K>({{2, 0}, {14, 0}, {2, 5}})), 2., 1, false);
    check_nodes_at_height("3-4-5 triangle", Polygon_with_holes_2(SST::make_polygon<K>({{0, 0}, {4, 0}, {0, 3}})), 1., 1, false);
    check_nodes_at_height("8-15-17 triangle", Polygon_with_holes_2(SST::make_polygon<K>({{0, 0}, {15, 0}, {0, 8}})), 3., 1, false);

    // Ridges (degenerate top contour): trapezoid with 3-4-5 sides, and rectangle
    check_nodes_at_height("3-4-5 trapezoid", Polygon_with_holes_2(SST::make_polygon<K>({{0, 0}, {10, 0}, {7, 4}, {3, 4}})), 2., 2, false);
    check_nodes_at_height("rectangle", Polygon_with_holes_2(SST::make_rectangle<K>(0, 0, 4, 2)), 1., 2, false);

    // A frame whose strips collapse at t = 0.5: the offset contour at that height has the nodes as corners;
    // then the corners of the left and right regions at t = 1 (degenerate top contour)
    {
      Polygon_with_holes_2 pwh(SST::make_rectangle<K>(0, 0, 8, 4));
      pwh.add_hole(SST::make_rectangle<K>(2, 1, 6, 3, false /*cw*/));
      check_nodes_at_height("frame (strips)", pwh, 0.5, 4, true);
      check_nodes_at_height("frame (corners)", pwh, 1., 4, false);
    }

    // House with a 3-4-5 roof: the bottom edge vanishes at t = 4.5 (the nodes are on the top contour,
    // which continues above), and the apex is at t = 5.25
    {
      const Polygon_with_holes_2 house(SST::make_polygon<K>({{0, 0}, {12, 0}, {12, 6}, {6, 14}, {0, 6}}));
      check_nodes_at_height("3-4-5 house (edge event)", house, 4.5, 2, true);
      check_nodes_at_height("3-4-5 house (apex)", house, 5.25, 1, false);
    }
  }

  // The height is not an event time: all points of the top contour are exactly at the height,
  // including those computed with the exact kernel fallback, for several heights that are not exactly
  // representable, or close to event times.
  static void test_top_contour_at_height()
  {
    const std::vector<Polygon_with_holes_2> inputs = {
      Polygon_with_holes_2(SST::make_polygon<K>({{0, 0}, {4, 3}, {1, 7}, {-3, 4}})),
      Polygon_with_holes_2(SST::make_polygon<K>({{0, 0}, {10, 0}, {7, 4}, {3, 4}})),
      Polygon_with_holes_2(SST::make_polygon<K>({{0, 0}, {4, 0}, {4, 2}, {2.5, 2}, {2, 2.5}, {2, 4}, {0, 4}})),
      Polygon_with_holes_2(SST::make_polygon<K>({{0.1, 0.2}, {3.3, 0.1}, {3.7, 2.9}, {1.9, 1.3}, {0.3, 3.1}}))
    };

    for(std::size_t i=0; i<inputs.size(); ++i)
    {
      for(double H : {0.1, 1. / 3., 0.7071067811865476, 1. - 1e-12, 1. + 1e-12})
      {
        for(const double sign : {1., -1.})
        {
          const std::string name = "top contour #" + std::to_string(i) + ", H = " + std::to_string(sign * H);
          Input in(inputs[i]);
          in.with_height(sign * H);
          Expectations e;
          e.rel_tol = 1e-6;
          const Mesh m = run(name, in, e);

          // - if some nodes are above the height, there is a top contour, exactly at the height
          // - unless a node is (legitimately) just below the height, there is nothing just below the height
          bool has_nodes_above = false, has_nodes_just_below = false;
          for(const auto& n : skeleton_nodes(inputs[i]))
          {
            if(n.second > EPECK::FT(H))
              has_nodes_above = true;
            else if(n.second >= EPECK::FT(H * (1 - 1e-9)))
              has_nodes_just_below = true;
          }

          std::size_t top_n = 0;
          for(vertex_descriptor v : vertices(m))
          {
            const double z = std::abs(CGAL::to_double(m.point(v).z()));
            if(z == H)
              ++top_n;
            if(!has_nodes_just_below)
              CGAL_SLS_CHECK_MSG(z == H || z < H * (1 - 1e-9), name << ": " << m.point(v) << " is just below the height");
          }

          if(has_nodes_above)
            CGAL_SLS_CHECK_MSG(top_n >= 3, name << ": only " << top_n << " points at the height");
          else
            CGAL_SLS_CHECK_MSG(top_n == 0, name << ": " << top_n << " points at the height, above the skeleton");
        }
      }
    }
  }

  // Snapping and clamping together: the top points on vertical walls are snapped in xy, and exactly at the height
  static void test_walls_at_height()
  {
    // all vertical except one edge with a very gentle slope: the top is a thin sliver near the sloped edge
    const Polygon_2 sq = SST::make_rectangle<K>(0, 0, 2, 2);
    for(double H : {0.5, 1. / 3., -0.7})
    {
      const std::string name = "three walls, H = " + std::to_string(H);
      const double h = std::abs(H);
      Input in(sq);
      in.with_weights(Weights{{FT(0), FT(0), FT(0), FT(1)}}).with_height(H); // the top edge (v2, v3) is sloped
      Expectations e;
      e.nv = 8;
      e.nf = 12;
      e.volume = 4 * h - h * h; // integral of 2(2 - z)
      e.top_area = 2 * (2 - h);
      e.rel_tol = 1e-6;
      const Mesh m = run(name, in, e);
      CGAL_SLS_CHECK_MSG(SST::has_vertex(m, Point_3(0, 0, H)), name << ": corner (0, 0)");
      CGAL_SLS_CHECK_MSG(SST::has_vertex(m, Point_3(2, 0, H)), name << ": corner (2, 0)");
      CGAL_SLS_CHECK_MSG(check_wall(name, m, P(0, 2), P(0, 0)) == 4, name << ": left wall");
      CGAL_SLS_CHECK_MSG(check_wall(name, m, P(2, 0), P(2, 2)) == 4, name << ": right wall");
    }
  }

  static void run_all()
  {
    // vertical slabs
    test_axis_aligned_wall();
    test_rotated_wall();
    test_adjacent_walls();
    test_non_adjacent_walls();
    test_hole_walls();
    test_outward_wall();
    test_prisms();
    // test_height_beyond_vertical_slabs();

    // height
    test_nodes_at_height();
    test_top_contour_at_height();

    // both
    test_walls_at_height();
  }
};

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  std::cout << "== EPICK ==" << std::endl;
  Snapping_tests<EPICK>::run_all();
  std::cout << "== EPECK ==" << std::endl;
  Snapping_tests<EPECK>::run_all();

  return SST::test_exit_code();
}
