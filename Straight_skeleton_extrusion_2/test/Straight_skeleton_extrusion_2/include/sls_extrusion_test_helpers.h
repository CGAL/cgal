#ifndef CGAL_SLS_EXTRUSION_TEST_HELPERS_H
#define CGAL_SLS_EXTRUSION_TEST_HELPERS_H

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

// #define CGAL_SLS_OUTPUT_FILES
// #define CGAL_SLS_DEBUG_DRAW

// #define CGAL_STRAIGHT_SKELETON_ENABLE_TRACE 2
// #define CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE
// #define CGAL_POLYGON_OFFSET_ENABLE_TRACE 4

void Straight_skeleton_external_trace ( std::string m )
{
  std::cout << std::setprecision(17) << m << std::endl;
}
void Straight_skeleton_traits_external_trace ( std::string m )
{
  std::cout << std::setprecision(17) << m << std::endl;
}

void error_handler ( char const* what, char const* expr, char const* file, int line, char const* msg )
{
  std::cerr << "CGAL error: " << what << " violation!" << std::endl
       << "Expr: " << expr << std::endl
       << "File: " << file << std::endl
       << "Line: " << line << std::endl;
  if ( msg != nullptr)
      std::cerr << "Explanation:" << msg << std::endl;

  std::exit(1);
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

#ifdef CGAL_SLS_DEBUG_DRAW
# include <CGAL/draw_surface_mesh.h>
# include <CGAL/draw_polyhedron.h>
# include <CGAL/draw_polygon_2.h>
# include <CGAL/draw_polygon_with_holes_2.h>
#endif

#include <CGAL/extrude_skeleton.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace CGAL {
namespace Straight_skeletons_2 {
namespace Extrusion_tests {

int& failure_count()
{
  static int n = 0;
  return n;
}

int test_exit_code()
{
  if(failure_count() != 0)
  {
    std::cerr << failure_count() << " check(s) FAILED" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Done" << std::endl;
  return EXIT_SUCCESS;
}

#define CGAL_SLS_CHECK(cond)                                                                       \
  do {                                                                                             \
    if(!(cond)) {                                                                                  \
      ++::CGAL::Straight_skeletons_2::Extrusion_tests::failure_count();                            \
      std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK FAILED: " #cond << std::endl;           \
    }                                                                                              \
  } while(0)

#define CGAL_SLS_CHECK_MSG(cond, msg)                                                              \
  do {                                                                                             \
    if(!(cond)) {                                                                                  \
      ++::CGAL::Straight_skeletons_2::Extrusion_tests::failure_count();                            \
      std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK FAILED: " #cond                         \
                << " [" << msg << "]" << std::endl;                                                \
    }                                                                                              \
  } while(0)

// Whether CGAL preconditions / assertions / postconditions are compiled in.
// Used to check invalid inputs.
#if defined(CGAL_NDEBUG) || (defined(NDEBUG) && !defined(CGAL_DEBUG))
constexpr bool has_cgal_preconditions = false;
constexpr bool has_cgal_assertions = false;
#else
# ifdef CGAL_NO_PRECONDITIONS
constexpr bool has_cgal_preconditions = false;
# else
constexpr bool has_cgal_preconditions = true;
# endif
# ifdef CGAL_NO_ASSERTIONS
constexpr bool has_cgal_assertions = false;
# else
constexpr bool has_cgal_assertions = true;
# endif
# ifdef CGAL_NO_POSTCONDITIONS
constexpr bool has_cgal_postconditions = false;
# else
constexpr bool has_cgal_postconditions = true;
# endif
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
// Input description

template <typename K>
struct Extrusion_input
{
  using FT = typename K::FT;
  using Weights = std::vector<std::vector<FT> >;

  CGAL::Polygon_with_holes_2<K> pwh;
  Weights weights; // empty: not passed
  Weights angles; // empty: not passed
  std::optional<FT> height; // nullopt: not passed

  Extrusion_input() = default;
  explicit Extrusion_input(const CGAL::Polygon_with_holes_2<K>& p) : pwh(p) { }
  explicit Extrusion_input(const CGAL::Polygon_2<K>& p) : pwh(p) { }

  Extrusion_input& with_weights(const Weights& w) { weights = w; return *this; }
  Extrusion_input& with_angles(const Weights& a) { angles = a; return *this; }
  Extrusion_input& with_height(const FT h) { height = h; return *this; }
  Extrusion_input& without_height() { height.reset(); return *this; }

  // Fill a weight range with the same value for every edge of every contour
  static Weights uniform(const CGAL::Polygon_with_holes_2<K>& p, const double v)
  {
    Weights w;
    w.emplace_back(p.outer_boundary().size(), FT(v));
    for(const auto& h : p.holes())
      w.emplace_back(h.size(), FT(v));
    return w;
  }
};

template <typename K>
CGAL::Polygon_2<K> make_polygon(std::initializer_list<std::array<double, 2> > pts)
{
  CGAL::Polygon_2<K> p;
  for(const auto& q : pts)
    p.push_back(typename K::Point_2(q[0], q[1]));
  return p;
}

// Axis-aligned rectangle; `ccw` selects the orientation (holes must be clockwise)
template <typename K>
CGAL::Polygon_2<K> make_rectangle(double x0, double y0, double x1, double y1, bool ccw = true)
{
  CGAL::Polygon_2<K> p = make_polygon<K>({{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}});
  if(!ccw)
    p.reverse_orientation();
  return p;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Running the function and classifying the outcome

enum class Outcome
{
  SUCCESS,                // returned 'true'
  RETURNED_FALSE,         // returned 'false'
  PRECONDITION_EXCEPTION, // CGAL::Precondition_exception
  ASSERTION_EXCEPTION,    // CGAL::Assertion_exception
  POSTCONDITION_EXCEPTION, // CGAL::Postcondition_exception
  OTHER_EXCEPTION
};

const char* to_string(const Outcome& o)
{
  switch(o)
  {
    case Outcome::SUCCESS:
      return "SUCCESS";
    case Outcome::RETURNED_FALSE:
      return "RETURNED_FALSE";
    case Outcome::PRECONDITION_EXCEPTION:
      return "PRECONDITION_EXCEPTION";
    case Outcome::ASSERTION_EXCEPTION:
      return "ASSERTION_EXCEPTION";
    case Outcome::POSTCONDITION_EXCEPTION:
      return "POSTCONDITION_EXCEPTION";
    default: return "OTHER_EXCEPTION";
  }
}

template <typename K, typename Mesh>
Outcome extrude(const Extrusion_input<K>& in,
                Mesh& out,
                const bool verbose = false)
{
  // the weights() and angles() named parameters take a non-const reference: use local copies
  auto weights = in.weights;
  auto angles = in.angles;

  try
  {
    bool res = false;
    const auto np = CGAL::parameters::verbose(verbose);
    if(!in.weights.empty())
    {
      if(in.height)
        res = CGAL::extrude_skeleton(in.pwh, out, np.weights(weights).maximum_height(*in.height));
      else
        res = CGAL::extrude_skeleton(in.pwh, out, np.weights(weights));
    }
    else if(!in.angles.empty())
    {
      if(in.height)
        res = CGAL::extrude_skeleton(in.pwh, out, np.angles(angles).maximum_height(*in.height));
      else
        res = CGAL::extrude_skeleton(in.pwh, out, np.angles(angles));
    }
    else
    {
      if(in.height)
        res = CGAL::extrude_skeleton(in.pwh, out, np.maximum_height(*in.height));
      else
        res = CGAL::extrude_skeleton(in.pwh, out, np);
    }
    return res ? Outcome::SUCCESS : Outcome::RETURNED_FALSE;
  }
  catch(const CGAL::Precondition_exception&) { return Outcome::PRECONDITION_EXCEPTION; }
  catch(const CGAL::Assertion_exception&) { return Outcome::ASSERTION_EXCEPTION; }
  catch(const CGAL::Postcondition_exception&) { return Outcome::POSTCONDITION_EXCEPTION; }
  catch(const std::exception&) { return Outcome::OTHER_EXCEPTION; }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Helpers

template <typename Mesh>
void flip_mesh(Mesh& m)
{
  using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;
  auto vpm = get(CGAL::vertex_point, m);
  using P = typename boost::property_traits<decltype(vpm)>::value_type;
  for (vertex_descriptor v : vertices(m))
  {
    const auto& p = get(vpm, v);
    put(vpm, v, P{p.x(), p.y(), -p.z()});
  }
  PMP::reverse_face_orientations(m);
}

// to choose a tolerance relative to the size of the input
double diag_xy(const CGAL::Bbox_3& b)
{
  return std::sqrt(CGAL::square(b.xmax() - b.xmin()) + CGAL::square(b.ymax() - b.ymin()));
}

template <typename K>
double polygon_area(const CGAL::Polygon_with_holes_2<K>& pwh)
{
  double a = CGAL::to_double(pwh.outer_boundary().area());
  for(const auto& h : pwh.holes())
    a += CGAL::to_double(h.area()); // the area is signed, hence negative for holes
  return a;
}

template <typename K, typename Mesh>
double clipped_volume(Mesh m /*copy*/, const double H)
{
  const double z = (H > 0) ? 1. : -1.;
  const typename K::Plane_3 pl { 0, 0, z, -H};
  PMP::clip(m, pl, CGAL::parameters::clip_volume(true));
  return CGAL::to_double(PMP::volume(m));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Checks

template <typename NT1, typename NT2>
bool approx_equal(const NT1& ia, const NT2& ib,
                  const double rel, const double abs_tol = 0.)
{
  const double a = CGAL::to_double(ia);
  const double b = CGAL::to_double(ib);
  const double diff = std::abs(a - b);
  return diff <= abs_tol || diff <= rel * (std::max)(std::abs(a), std::abs(b));
}

template <typename Mesh, typename P>
bool has_vertex(const Mesh& m, const P& p)
{
  using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;
  for(vertex_descriptor v : vertices(m))
    if(m.point(v) == p)
      return true;
  return false;
}

template <typename Mesh, typename P>
bool has_vertex_approx(const Mesh& m, const P& p, double tol = 1e-9)
{
  using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;
  for(vertex_descriptor v : vertices(m))
  {
    const P& q = m.point(v);
    if(CGAL::abs(q.x() - p.x()) <= tol &&
       CGAL::abs(q.y() - p.y()) <= tol &&
       CGAL::abs(q.z() - p.z()) <= tol)
      return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

struct Expectations
{
  std::optional<std::size_t> nv, nf;   // expected #vertices / #faces (triangles)
  std::optional<double> volume;        // expected |volume|
  std::optional<double> top_area;      // expected area of the top horizontal set of faces
  std::optional<double> zmax_abs;      // expected max |z|
  std::optional<int> genus;            // expected genus of the closed output
  bool check_self_intersections = true;
  double rel_tol = 1e-7;               // relative tolerance for geometric comparisons
};

template <typename K, typename Mesh>
void check_extrusion(const std::string& name,
                     const Extrusion_input<K>& in,
                     const Mesh& m,
                     const Expectations& e = Expectations())
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  using FT = typename K::FT;
  using Point_3 = typename K::Point_3;

  using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;
  using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;

  const int failures_n = failure_count();

  auto vpm = get(CGAL::vertex_point, m);

  const bool capped = in.height.has_value();
  const double dir = (in.height && CGAL::is_negative(*in.height)) ? -1. : 1.;

  // 1. Mesh validity
  CGAL_SLS_CHECK_MSG(num_faces(m) > 0, name);
  CGAL_SLS_CHECK_MSG(CGAL::is_valid_polygon_mesh(m), name);
  CGAL_SLS_CHECK_MSG(CGAL::is_closed(m), name);
  CGAL_SLS_CHECK_MSG(CGAL::is_triangle_mesh(m), name);
  if(!CGAL::is_closed(m) || !CGAL::is_triangle_mesh(m))
    return;

  CGAL_SLS_CHECK_MSG(PMP::is_outward_oriented(m), name);

  // 2. Combinatorial info
  if(e.nv)
  {
    CGAL_SLS_CHECK_MSG(num_vertices(m) == *e.nv,
                       name << ": #V = " << num_vertices(m) << " expected " << *e.nv);
  }

  if(e.nf)
  {
    CGAL_SLS_CHECK_MSG(num_faces(m) == *e.nf,
                       name << ": #F = " << num_faces(m) << " expected " << *e.nf);
  }

  if(e.genus)
  {
    const long chi = long(num_vertices(m)) - long(num_edges(m)) + long(num_faces(m));
    CGAL_SLS_CHECK_MSG(chi == 2 - 2 * (*e.genus), name << ": Euler characteristic " << chi);
  }

  // 3. Geometric embedding
  if (CGAL::Straight_skeleton_extrusion::internal::is_exact_kernel<K>)
  {
    std::size_t degenerate_faces = 0;
    for(face_descriptor f : faces(m))
      if(PMP::is_degenerate_triangle_face(f, m))
        ++degenerate_faces;
    CGAL_SLS_CHECK_MSG(degenerate_faces == 0, name << ": " << degenerate_faces << " degenerate faces");
  }

  if(e.check_self_intersections)
  {
    if (CGAL::Straight_skeleton_extrusion::internal::is_exact_kernel<K>)
      CGAL_SLS_CHECK_MSG(!PMP::does_self_intersect(m), name);
    else if (PMP::does_self_intersect(m))
      std::cerr << "W: Self-intersections with kernel " << typeid(K).name() << std::endl;
  }

  // 4. Height
  const CGAL::Bbox_3 bb = PMP::bbox(m);

  for(vertex_descriptor v : vertices(m))
  {
    const FT z = get(vpm, v).z();
    if(dir > 0)
      CGAL_SLS_CHECK_MSG(!CGAL::is_negative(z), name << ": z = " << z);
    else
      CGAL_SLS_CHECK_MSG(!CGAL::is_positive(z), name << ": z = " << z);

    if(capped) // exact comparison: heights are clamped to the maximum height
      CGAL_SLS_CHECK_MSG(CGAL::abs(z) <= CGAL::abs(*in.height),
                         name << ": z = " << z << " above max height");
  }

  if(e.zmax_abs)
  {
    const double zm = (dir > 0) ? bb.zmax() : -bb.zmin();
    CGAL_SLS_CHECK_MSG(approx_equal(zm, *e.zmax_abs, e.rel_tol, e.rel_tol * (*e.zmax_abs)),
                       name << ": max |z| = " << zm << " expected " << *e.zmax_abs);
  }

  // 4bis.
  double top_area = 0.;
  for(face_descriptor f : faces(m))
  {
    std::array<Point_3, 3> t;
    int i = 0;
    for(vertex_descriptor v : CGAL::vertices_around_face(halfedge(f, m), m))
    {
      const auto& p = get(vpm, v);
      t[i++] = p;
    }
    const bool is_top = capped && t[0][2] == *in.height && t[1][2] == *in.height && t[2][2] == *in.height;
    if(is_top) {
      // don't use PMP::face_area() because it would bring a sqrt, and we know we are
      // in a z-constant face
      top_area += 0.5 * CGAL::to_double((t[1].x() - t[0].x()) * (t[2].y() - t[0].y()) - (t[1].y() - t[0].y()) * (t[2].x() - t[0].x()));
    }
  }

  const double scale = diag_xy(bb);
  const double atol = e.rel_tol * scale * scale;
  if(e.top_area)
  {
    CGAL_SLS_CHECK_MSG(approx_equal(dir * top_area, *e.top_area, e.rel_tol, atol),
                       name << ": top cap area " << top_area << " expected " << *e.top_area);
  }
  else if(capped)
  {
    CGAL_SLS_CHECK_MSG(!CGAL::is_negative(FT(dir * top_area)),
                       name << ": top cap is facing the wrong way");
  }

  // 5. Volume
  const FT vol = PMP::volume(m);
  CGAL_SLS_CHECK_MSG(vol > 0, name << ": volume " << vol);
  if(e.volume)
  {
    CGAL_SLS_CHECK_MSG(approx_equal(vol, *e.volume, e.rel_tol, e.rel_tol * scale * scale * scale),
                       name << ": volume " << vol << " expected " << *e.volume);

  }

  if(failure_count() != failures_n)
    std::cerr << "  ... in test '" << name << "'" << std::endl;
}

} // namespace Extrusion_tests
} // namespace Straight_skeletons_2

namespace Polygon_mesh_processing {

// preview of https://github.com/CGAL/cgal/pull/9828
template <typename PolygonMesh1,
          typename PolygonMesh2,
          typename NamedParameters1 = parameters::Default_named_parameters,
          typename NamedParameters2 = parameters::Default_named_parameters>
bool are_equal(const PolygonMesh1& m1,
               const PolygonMesh2& m2,
               const NamedParameters1& np1 = parameters::default_values(),
               const NamedParameters2& np2 = parameters::default_values())
{
  if (vertices(m1).size() != vertices(m2).size() ||
      edges(m1).size()    != edges(m2).size()    ||
      faces(m1).size()    != faces(m2).size())
    return false;

  auto throw_at_output = boost::make_function_output_iterator(CGAL::internal::Throw_at_output{});
  try
  {
    CGAL::Polygon_mesh_processing::match_faces(m1, m2, CGAL::Emptyset_iterator(), throw_at_output, throw_at_output, np1, np2);
  }
  catch (const CGAL::internal::Throw_at_output_exception&)
  {
    return false;
  }

  return true;
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_SLS_EXTRUSION_TEST_HELPERS_H
