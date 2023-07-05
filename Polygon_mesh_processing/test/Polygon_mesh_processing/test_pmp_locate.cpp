#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

// Graphs
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Regular_triangulation_2.h>
#include <CGAL/boost/graph/properties_Regular_triangulation_2.h>

#include <CGAL/Polygon_mesh_processing/locate.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Dimension.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Origin.h>
#include <CGAL/property_map.h>
#include <CGAL/Random.h>
#include <CGAL/Unique_hash_map.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <optional>

#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel    EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel      EPECK;

typedef CGAL::Simple_cartesian<typename CGAL::internal::Exact_field_selector<double>::Type> Exact_kernel;

template<typename AABB_tree>
typename CGAL::Kernel_traits<typename AABB_tree::AABB_traits::Point_3>::type::Ray_2
random_2D_ray(const AABB_tree& aabb_tree, CGAL::Random& rnd)
{
  typedef typename AABB_tree::AABB_traits::Point_3           Point_3;
  typedef typename CGAL::Kernel_traits<Point_3>::type        Kernel;
  typedef typename Kernel::FT                                FT;
  typedef typename Kernel::Point_2                           Point_2;
  typedef typename Kernel::Ray_2                             Ray_2;

  const CGAL::Bbox_3& bbox = aabb_tree.bbox();

  FT px = (bbox.xmin() == bbox.xmax()) ? bbox.xmin() : rnd.get_double(bbox.xmin(), bbox.xmax());
  FT py = (bbox.ymin() == bbox.ymax()) ? bbox.ymin() : rnd.get_double(bbox.ymin(), bbox.ymax());

  FT qx = (bbox.xmin() == bbox.xmax()) ? bbox.xmin() : rnd.get_double(bbox.xmin(), bbox.xmax());
  FT qy = (bbox.ymin() == bbox.ymax()) ? bbox.ymin() : rnd.get_double(bbox.ymin(), bbox.ymax());

  return Ray_2(Point_2(px, py), Point_2(qx, qy));
}

template<typename AABB_tree>
typename CGAL::Kernel_traits<typename AABB_tree::AABB_traits::Point_3>::type::Ray_3
random_3D_ray(const AABB_tree& aabb_tree, CGAL::Random& rnd)
{
  typedef typename AABB_tree::AABB_traits::Point_3           Point_3;
  typedef typename CGAL::Kernel_traits<Point_3>::type        Kernel;
  typedef typename Kernel::FT                                FT;
  typedef typename Kernel::Ray_3                             Ray_3;

  const CGAL::Bbox_3& bbox = aabb_tree.bbox();

  FT px = (bbox.xmin() == bbox.xmax()) ? bbox.xmin() : rnd.get_double(bbox.xmin(), bbox.xmax());
  FT py = (bbox.ymin() == bbox.ymax()) ? bbox.ymin() : rnd.get_double(bbox.ymin(), bbox.ymax());
  FT pz = (bbox.zmin() == bbox.zmax()) ? bbox.zmin() : rnd.get_double(bbox.zmin(), bbox.zmax());

  FT qx = (bbox.xmin() == bbox.xmax()) ? bbox.xmin() : rnd.get_double(bbox.xmin(), bbox.xmax());
  FT qy = (bbox.ymin() == bbox.ymax()) ? bbox.ymin() : rnd.get_double(bbox.ymin(), bbox.ymax());
  FT qz = (bbox.zmin() == bbox.zmax()) ? bbox.zmin() : rnd.get_double(bbox.zmin(), bbox.zmax());

  return Ray_3(Point_3(px, py, pz), Point_3(qx, qy, qz));
}

template<typename FT>
bool is_equal(const FT& a, const FT& b)
{
  if(std::is_floating_point<FT>::value)
    return (CGAL::abs(a - b) <= 1e-7); // numeric_limits' epsilon is too restrictive...
  else
    return (a == b);
}

template<typename K, typename G>
void test_snappers(const G& g)
{
  std::cout << "  test snappers..." << std::endl;

  typedef typename K::FT                                              FT;

  PMP::Barycentric_coordinates<FT> coords = CGAL::make_array(FT(1e-6), FT(0.9999999999999999999), FT(1e-7));
  PMP::Face_location<G, FT> loc = std::make_pair(*(faces(g).first), coords);

  // ---------------------------------------------------------------------------
  PMP::internal::snap_coordinates_to_border<FT>(coords); // uses numeric_limits' epsilon()
  assert(coords[0] == FT(1e-6) && coords[1] == FT(1) && coords[2] == FT(1e-7));

  PMP::internal::snap_coordinates_to_border(coords, FT(1e-5));
  assert(coords[0] == FT(0) && coords[1] == FT(1) && coords[2] == FT(0));

  // ---------------------------------------------------------------------------
  PMP::internal::snap_location_to_border<FT>(loc, g); // uses numeric_limits' epsilon()
  assert(!PMP::is_on_face_border(loc, g));

  PMP::internal::snap_location_to_border(loc, g, FT(1e-7));
  assert(PMP::is_on_face_border(loc, g));
}

template <typename K, int>
struct Point_to_bare_point
{
  typedef typename K::Point_2 type;
};

template <typename K>
struct Point_to_bare_point<K, 3>
{
  typedef typename K::Point_3 type;
};

template<typename K, typename G, typename VPM>
void test_constructions(const G& g,
                        const VPM vpm,
                        CGAL::Random& rnd)
{
  std::cout << "  test constructions..." << std::endl;

  typedef typename boost::graph_traits<G>::vertex_descriptor                 vertex_descriptor;
  typedef typename boost::graph_traits<G>::halfedge_descriptor               halfedge_descriptor;
  typedef typename boost::graph_traits<G>::face_descriptor                   face_descriptor;
  typedef typename PMP::descriptor_variant<G>                                descriptor_variant;

  typedef typename boost::property_traits<VPM>::value_type                        Point;
  typedef typename boost::property_traits<VPM>::reference                         Point_reference;
  typedef typename K::FT                                                          FT;
  typedef typename Point_to_bare_point<K, Point::Ambient_dimension::value>::type  Bare_point;

  typedef typename PMP::Barycentric_coordinates<FT>                          Barycentric_coordinates;
  typedef typename PMP::Face_location<G, FT>                                 Face_location;

  face_descriptor f = CGAL::internal::random_face_in_mesh(g, rnd);
  halfedge_descriptor h = halfedge(f, g);
  vertex_descriptor v = source(h, g);

  Point_reference p = get(vpm, v);
  Point_reference q = get(vpm, target(h, g));
  Point_reference r = get(vpm, target(next(h, g), g));

  const Bare_point bp(p);
  const Bare_point bq(q);
  const Bare_point br(r);

  Barycentric_coordinates bar;
  Face_location loc;
  loc.first = f;

  // ---------------------------------------------------------------------------
  bar = PMP::barycentric_coordinates(p, q, r, p, K());
  assert(is_equal(bar[0], FT(1)) && is_equal(bar[1], FT(0)) && is_equal(bar[2], FT(0)));

  bar = PMP::barycentric_coordinates(p, q, r, q, K());
  assert(is_equal(bar[0], FT(0)) && is_equal(bar[1], FT(1)) && is_equal(bar[2], FT(0)));

  bar = PMP::barycentric_coordinates(p, q, r, r, K());
  assert(is_equal(bar[0], FT(0)) && is_equal(bar[1], FT(0)) && is_equal(bar[2], FT(1)));

  Point mp = Point(CGAL::midpoint(bp, bq));
  bar = PMP::barycentric_coordinates(p, q, r, mp);
  assert(is_equal(bar[0], FT(0.5)) && is_equal(bar[1], FT(0.5)) && is_equal(bar[2], FT(0)));

  int n = 100;
  while(n --> 0) // :)
  {
    const FT a = rnd.get_double(-1., 1.);
    const FT b = rnd.get_double(-1., 1.);
    const FT c = 1. - a - b;

    // Point to location and inversely
    Bare_point barycentric_pt = CGAL::barycenter(bp, a, bq, b, br, c);
    bar = PMP::barycentric_coordinates(p, q, r, Point(barycentric_pt));
    assert(is_equal(bar[0], a) && is_equal(bar[1], b) && is_equal(bar[2], c));

    loc.second = bar;
    const Bare_point barycentric_pt_2 =
      Bare_point(PMP::construct_point(loc, g,
                                      CGAL::parameters::vertex_point_map(vpm)
                                                       .geom_traits(K())));

    const FT sq_dist = CGAL::squared_distance(barycentric_pt, barycentric_pt_2);
    assert(is_equal(sq_dist, FT(0)));
  }

  // ---------------------------------------------------------------------------
  loc = std::make_pair(f, CGAL::make_array(FT(0.3), FT(0.4), FT(0.3)));
  descriptor_variant dv = PMP::get_descriptor_from_location(loc, g);
  const face_descriptor* fd = std::get_if<face_descriptor>(&dv);
  assert(fd);

  loc = std::make_pair(f, CGAL::make_array(FT(0.5), FT(0.5), FT(0)));
  dv = PMP::get_descriptor_from_location(loc, g);
  const halfedge_descriptor* hd = std::get_if<halfedge_descriptor>(&dv);
  assert(hd);

  loc = std::make_pair(f, CGAL::make_array(FT(1), FT(0), FT(0)));
  assert(PMP::is_on_vertex(loc, source(halfedge(f, g), g), g));

  dv = PMP::get_descriptor_from_location(loc, g);
  if(const vertex_descriptor* v = std::get_if<vertex_descriptor>(&dv)) { } else { assert(false); }

  // ---------------------------------------------------------------------------
  // just to check the API
  PMP::construct_point(loc, g);
  PMP::construct_point(loc, g, CGAL::parameters::default_values());
}

template<typename K, typename G>
void test_random_entities(const G& g, CGAL::Random& rnd)
{
  std::cout << "  test random entities..." << std::endl;

  typedef typename boost::graph_traits<G>::halfedge_descriptor               halfedge_descriptor;
  typedef typename boost::graph_traits<G>::face_descriptor                   face_descriptor;

  typedef typename K::FT                                                     FT;
  typedef typename PMP::Face_location<G, FT>                                 Face_location;

  // ---------------------------------------------------------------------------
  Face_location loc;

  halfedge_descriptor h = CGAL::internal::random_halfedge_in_mesh(g, rnd);
  if(is_border(h, g))
    h = opposite(h, g);
  face_descriptor f = CGAL::internal::random_face_in_mesh(g, rnd);

  int nn = 100;
  while(nn --> 0) // the infamous 'go to zero' operator
  {
    loc = PMP::random_location_on_mesh<FT>(g, rnd);
    assert(loc.first != boost::graph_traits<G>::null_face());
    assert(loc.second[0] >= FT(0) && loc.second[0] <= FT(1) &&
      loc.second[1] >= FT(0) && loc.second[1] <= FT(1) &&
      loc.second[2] >= FT(0) && loc.second[2] <= FT(1));

    loc = PMP::random_location_on_face<FT>(f, g, rnd);
    assert(loc.first == f);
    assert(loc.second[0] >= FT(0) && loc.second[0] <= FT(1) &&
      loc.second[1] >= FT(0) && loc.second[1] <= FT(1) &&
      loc.second[2] >= FT(0) && loc.second[2] <= FT(1));

    loc = PMP::random_location_on_halfedge<FT>(h, g, rnd);
    assert(loc.first == face(h, g));
    assert(loc.second[0] >= FT(0) && loc.second[0] <= FT(1) &&
      loc.second[1] >= FT(0) && loc.second[1] <= FT(1) &&
      loc.second[2] >= FT(0) && loc.second[2] <= FT(1));
    int h_id = CGAL::halfedge_index_in_face(h, g);

    assert(loc.second[(h_id + 2) % 3] == FT(0));
  }
}

template<typename K, typename G>
void test_helpers(const G& g, CGAL::Random& rnd)
{
  std::cout << "  test helpers..." << std::endl;

  typedef typename boost::graph_traits<G>::vertex_descriptor                 vertex_descriptor;
  typedef typename boost::graph_traits<G>::halfedge_descriptor               halfedge_descriptor;
  typedef typename boost::graph_traits<G>::face_descriptor                   face_descriptor;

  typedef typename K::FT                                                     FT;
  typedef typename PMP::Face_location<G, FT>                                 Face_location;

  face_descriptor f = CGAL::internal::random_face_in_mesh(g, rnd);
  halfedge_descriptor h = halfedge(f, g);
  vertex_descriptor v = source(h, g);

  // ---------------------------------------------------------------------------
  // Local index
  int pos = CGAL::vertex_index_in_face(v, f, g);
  assert(pos == 0);
  pos = CGAL::vertex_index_in_face(target(h, g), f, g);
  assert(pos == 1);
  pos = CGAL::vertex_index_in_face(target(next(h, g), g), f, g);
  assert(pos == 2);

  pos = CGAL::halfedge_index_in_face(h, g);
  assert(pos == 0);
  pos = CGAL::halfedge_index_in_face(next(h, g), g);
  assert(pos == 1);
  pos = CGAL::halfedge_index_in_face(prev(h, g), g);
  assert(pos == 2);

  // ---------------------------------------------------------------------------
  // Incident faces
  Face_location loc = PMP::random_location_on_face<FT>(f, g, rnd);
  std::set<face_descriptor> s;
  PMP::internal::incident_faces(loc, g, std::inserter(s, s.begin()));
  assert(PMP::is_on_face_border(loc, g) || s.size() == 1);

  loc = PMP::random_location_on_halfedge<FT>(h, g, rnd);
  std::vector<face_descriptor> vec;
  PMP::internal::incident_faces(loc, g, std::back_inserter(vec));
  assert(PMP::is_on_vertex(loc, source(h, g), g) ||
    PMP::is_on_vertex(loc, target(h, g), g) ||
    vec.size() == 2);
}

template<typename K, typename G>
void test_predicates(const G& g, CGAL::Random& rnd)
{
  std::cout << "  test predicates..." << std::endl;

  typedef typename K::FT                                                     FT;

  typedef typename boost::graph_traits<G>::vertex_descriptor                 vertex_descriptor;
  typedef typename boost::graph_traits<G>::halfedge_descriptor               halfedge_descriptor;
  typedef typename boost::graph_traits<G>::face_descriptor                   face_descriptor;

  typedef typename PMP::Face_location<G, FT>                                 Face_location;

  face_descriptor f = CGAL::internal::random_face_in_mesh(g, rnd);
  halfedge_descriptor h = halfedge(f, g);
  vertex_descriptor v = source(h, g);

  // ---------------------------------------------------------------------------
  Face_location loc(f, CGAL::make_array(FT(1), FT(0), FT(0)));
  assert(PMP::is_on_vertex<FT>(loc, v, g));
  loc = Face_location(f, CGAL::make_array(FT(0), FT(1), FT(0)));
  assert(PMP::is_on_vertex<FT>(loc, target(h, g), g));
  loc = Face_location(f, CGAL::make_array(FT(0), FT(0), FT(1)));
  assert(PMP::is_on_vertex<FT>(loc, target(next(h, g), g), g));
  loc = Face_location(f, CGAL::make_array(FT(-1), FT(1), FT(1)));
  assert(!PMP::is_on_vertex<FT>(loc, target(next(h, g), g), g));

  // ---------------------------------------------------------------------------
  loc = Face_location(f, CGAL::make_array(FT(0.5), FT(0.5), FT(0)));
  assert(PMP::is_on_halfedge<FT>(loc, h, g));
  loc = Face_location(f, CGAL::make_array(FT(0), FT(0.5), FT(0.5)));
  assert(PMP::is_on_halfedge<FT>(loc, next(h, g), g));
  loc = Face_location(f, CGAL::make_array(FT(-0.5), FT(1.5), FT(0)));
  assert(!PMP::is_on_halfedge<FT>(loc, h, g));
  loc = Face_location(f, CGAL::make_array(FT(0.1), FT(-0.6), FT(1.5)));
  assert(!PMP::is_on_halfedge<FT>(loc, h, g));

  // ---------------------------------------------------------------------------
  loc = Face_location(f, CGAL::make_array(FT(0.3), FT(0.3), FT(0.4)));
  assert(PMP::is_in_face<FT>(loc, g));
  loc = Face_location(f, CGAL::make_array(FT(0), FT(0), FT(1)));
  assert(PMP::is_in_face<FT>(loc, g));
  loc = Face_location(f, CGAL::make_array(FT(0), FT(2), FT(-1)));
  assert(!PMP::is_in_face<FT>(loc, g));

  // ---------------------------------------------------------------------------
  loc = Face_location(f, CGAL::make_array(FT(0.3), FT(0.3), FT(0.4)));
  assert(!PMP::is_on_face_border<FT>(loc, g));
  loc = Face_location(f, CGAL::make_array(FT(0), FT(0.6), FT(0.4)));
  assert(PMP::is_on_face_border<FT>(loc, g));
  loc = Face_location(f, CGAL::make_array(FT(0), FT(0), FT(1)));
  assert(PMP::is_on_face_border<FT>(loc, g));
  loc = Face_location(f, CGAL::make_array(FT(-0.2), FT(0), FT(1.2)));
  assert(!PMP::is_on_face_border(loc, g));

  // ---------------------------------------------------------------------------
  int max = 1000, counter = 0;
  typename boost::graph_traits<G>::halfedge_iterator hit, hend;
  boost::tie(hit, hend) = halfedges(g);
  for(; hit!=hend; ++hit)
  {
    const halfedge_descriptor h = *hit;
    if(face(h, g) == boost::graph_traits<G>::null_face())
      continue;

    const int id_of_h = CGAL::halfedge_index_in_face(h, g);
    const face_descriptor f = face(h, g);
    loc.first = f;

    loc.second[id_of_h] = FT(1);
    loc.second[(id_of_h+1)%3] = FT(0);
    loc.second[(id_of_h+2)%3] = FT(0);
    std::optional<halfedge_descriptor> opt_hd = CGAL::is_border(source(h, g), g);
    assert(PMP::is_on_mesh_border<FT>(loc, g) == (opt_hd != std::nullopt));

    loc.second[id_of_h] = FT(0.5);
    loc.second[(id_of_h+1)%3] = FT(0.5);
    assert(PMP::is_on_mesh_border<FT>(loc, g) == CGAL::is_border(edge(h, g), g));

    // Even if the point does lie on the border of the mesh, 'false' is returned because
    // another face descriptor should be used.
    loc.second[id_of_h] = -0.5;
    loc.second[(id_of_h+1)%3] = 1.5;
    assert(!PMP::is_on_mesh_border<FT>(loc, g));

    if(++counter > max)
      break;
  }
}

template<typename K, typename G, typename VPM>
void test_locate_in_face(const G& g,
                         const VPM vpm,
                         CGAL::Random& rnd)
{
  std::cout << "  test locate_in_face()..." << std::endl;

  typedef typename boost::property_traits<VPM>::reference                    Point_reference;
  typedef typename K::FT                                                     FT;

  typedef typename boost::graph_traits<G>::vertex_descriptor                 vertex_descriptor;
  typedef typename boost::graph_traits<G>::halfedge_descriptor               halfedge_descriptor;
  typedef typename boost::graph_traits<G>::face_descriptor                   face_descriptor;

  typedef typename PMP::Face_location<G, FT>                                 Face_location;

  const face_descriptor f = CGAL::internal::random_face_in_mesh(g, rnd);
  const halfedge_descriptor h = halfedge(f, g);
  const vertex_descriptor v = target(h, g);

  Face_location loc;
  FT a = 0.1;
  Point_reference p = get(vpm, v);

  loc = PMP::locate_vertex<FT>(v, g);

  assert(is_equal(loc.second[CGAL::vertex_index_in_face(v, loc.first, g)], FT(1)));
  assert(is_equal(loc.second[(CGAL::vertex_index_in_face(v, loc.first, g) + 1) % 3], FT(0)));
  assert(is_equal(loc.second[(CGAL::vertex_index_in_face(v, loc.first, g) + 2) % 3], FT(0)));

  loc = PMP::locate_vertex<FT>(v, f, g);
  assert(loc.first == f);
  assert(is_equal(loc.second[0], FT(0)) && is_equal(loc.second[1], FT(1)) && is_equal(loc.second[2], FT(0)));

  loc = PMP::locate_on_halfedge<FT>(h, a, g);
  const int h_id = CGAL::halfedge_index_in_face(h, g);
  assert(loc.first == f && is_equal(loc.second[(h_id + 2) % 3], FT(0)));

  loc = PMP::locate_in_face(p, f, g, CGAL::parameters::vertex_point_map(vpm).geom_traits(K()));
  int v_id = CGAL::vertex_index_in_face(v, f, g);
  assert(loc.first == f && is_equal(loc.second[v_id], FT(1)));

  // Internal vertex point pmap
  typedef typename boost::property_map_value<G, CGAL::vertex_point_t>::type     Point;

  Point p2 = get(CGAL::vertex_point, g, v);
  PMP::locate_in_face(p2, f, g);
  assert(loc.first == f && is_equal(loc.second[v_id], FT(1)));

  // ---------------------------------------------------------------------------
  loc.second[0] = FT(0.2);
  loc.second[1] = FT(0.8);
  loc.second[2] = FT(0);

  halfedge_descriptor neigh_hd = opposite(halfedge(f, g), g);
  face_descriptor neigh_f = face(neigh_hd, g);

  // Want to check good correspondence seen from one side and the other. If unfortunately
  // we have selected a border face, can't do anything!
  if(neigh_f != boost::graph_traits<G>::null_face())
  {
    int neigh_hd_id = CGAL::halfedge_index_in_face(neigh_hd, g);
    Face_location neigh_loc;
    neigh_loc.first = neigh_f;
    neigh_loc.second[neigh_hd_id] = FT(0.3);
    neigh_loc.second[(neigh_hd_id+1)%3] = FT(0.7);
    neigh_loc.second[(neigh_hd_id+2)%3] = FT(0);

    PMP::locate_in_adjacent_face(loc, neigh_f, g);
    assert(PMP::locate_in_common_face<FT>(loc, neigh_loc, g));

    if (std::is_same<K, EPECK>()) {
      assert(PMP::locate_in_common_face<FT>(loc, p, neigh_loc, g, CGAL::parameters::vertex_point_map(vpm).geom_traits(K())));
    }
    assert(PMP::locate_in_common_face<FT>(loc, p, neigh_loc, g, CGAL::parameters::vertex_point_map(vpm).geom_traits(K()), 1e-7));
  }
}

template <typename K, typename VPM,
          int dim = CGAL::Ambient_dimension<typename boost::property_traits<VPM>::value_type>::value>
struct Locate_with_AABB_tree_Tester // 2D case
{
  template <typename G>
  void test(const G& g, const VPM vpm, CGAL::Random& rnd) const
  {
    std::cout << "  test locate_with_AABB_tree (2D)..." << std::endl;

    typedef typename boost::property_traits<VPM>::reference                    Point_reference;

    typedef typename K::FT                                                     FT;
    typedef typename K::Ray_2                                                  Ray_2;
    typedef typename K::Ray_3                                                  Ray_3;
    typedef typename K::Point_3                                                Point_3;

    typedef typename boost::graph_traits<G>::vertex_descriptor                 vertex_descriptor;
    typedef typename boost::graph_traits<G>::halfedge_descriptor               halfedge_descriptor;
    typedef typename boost::graph_traits<G>::face_descriptor                   face_descriptor;

    typedef typename PMP::Face_location<G, FT>                                 Face_location;

    face_descriptor f = CGAL::internal::random_face_in_mesh(g, rnd);
    halfedge_descriptor h = halfedge(f, g);
    vertex_descriptor v = target(h, g);

    // ---------------------------------------------------------------------------
    typedef typename boost::property_traits<VPM>::value_type                   Intrinsic_point;
    typedef PMP::internal::Point_to_Point_3<G, Intrinsic_point>                Intrinsic_point_to_Point_3;
    typedef PMP::internal::Point_to_Point_3_VPM<G, VPM>                        WrappedVPM;
    typedef CGAL::AABB_face_graph_triangle_primitive<G, WrappedVPM>            AABB_face_graph_primitive;
    typedef CGAL::AABB_traits<K, AABB_face_graph_primitive>                    AABB_face_graph_traits;

    static_assert(std::is_same<typename AABB_face_graph_traits::Point_3, Point_3>::value);

    Intrinsic_point_to_Point_3 to_p3;

    CGAL::AABB_tree<AABB_face_graph_traits> tree_a;
    Point_reference p_a = get(vpm, v);
    const Point_3& p3_a = to_p3(p_a);

    CGAL::AABB_tree<AABB_face_graph_traits> tree_b;
    WrappedVPM vpm_b(vpm);
    // ---------------------------------------------------------------------------

    PMP::build_AABB_tree(g, tree_a, CGAL::parameters::vertex_point_map(vpm));
    PMP::build_AABB_tree(g, tree_b, CGAL::parameters::vertex_point_map(vpm_b));
    assert(tree_b.size() == num_faces(g));

    Face_location loc = PMP::locate_with_AABB_tree(p_a, tree_a, g, CGAL::parameters::vertex_point_map(vpm));

    // sanitize otherwise some test platforms fail
    PMP::internal::snap_location_to_border(loc, g, FT(1e-7));
    assert(PMP::is_on_vertex(loc, v, g)); // might fail due to precision issues...
    assert(is_equal(loc.second[CGAL::vertex_index_in_face(v, loc.first, g)], FT(1)));
    assert(is_equal(loc.second[(CGAL::vertex_index_in_face(v, loc.first, g) + 1) % 3], FT(0)));
    assert(is_equal(loc.second[(CGAL::vertex_index_in_face(v, loc.first, g) + 2) % 3], FT(0)));
    assert(is_equal(CGAL::squared_distance(to_p3(
      PMP::construct_point<FT>(loc, g, CGAL::parameters::vertex_point_map(vpm))), p3_a), FT(0)));

    loc = PMP::locate_with_AABB_tree(p_a, tree_a, g, CGAL::parameters::vertex_point_map(vpm));
    assert(is_equal(CGAL::squared_distance(to_p3(
      PMP::construct_point<FT>(loc, g, CGAL::parameters::vertex_point_map(vpm))), p3_a), FT(0)));

    // ---------------------------------------------------------------------------
    loc = PMP::locate(p_a, g, CGAL::parameters::vertex_point_map(vpm));
    assert(is_equal(CGAL::squared_distance(to_p3(
      PMP::construct_point(loc, g, CGAL::parameters::vertex_point_map(vpm))), p3_a), FT(0)));

    if (std::is_same<K, EPECK>()) {
      assert(PMP::is_in_face(loc, g));
    }

    loc = PMP::locate_with_AABB_tree(CGAL::ORIGIN, tree_b, g, CGAL::parameters::vertex_point_map(vpm_b));
    if (std::is_same<K, EPECK>()) {
      assert(PMP::is_in_face(loc, g));
    }

    loc = PMP::locate(CGAL::ORIGIN, g, CGAL::parameters::vertex_point_map(vpm_b));
    if (std::is_same<K, EPECK>()) {
      assert(PMP::is_in_face(loc, g));
    }

    // ---------------------------------------------------------------------------
    Ray_2 r2 = random_2D_ray<CGAL::AABB_tree<AABB_face_graph_traits> >(tree_a, rnd);
    loc = PMP::locate_with_AABB_tree(r2, tree_a, g, CGAL::parameters::vertex_point_map(vpm));
    if(loc.first != boost::graph_traits<G>::null_face() && std::is_same<K, EPECK>())
      assert(PMP::is_in_face(loc, g));

    Ray_3 r3 = random_3D_ray<CGAL::AABB_tree<AABB_face_graph_traits> >(tree_b, rnd);
    loc = PMP::locate_with_AABB_tree(r3, tree_b, g, CGAL::parameters::vertex_point_map(vpm_b)
                                                                     .geom_traits(K()));
  }
};

template <typename K>
struct My_3D_Point
{
  typedef typename K::FT FT;

  typedef K R; // so that we can use Kernel_traits
  typedef CGAL::Dimension_tag<3>  Ambient_dimension;
  typedef CGAL::Dimension_tag<0>  Feature_dimension;

  My_3D_Point() { }
  My_3D_Point(const CGAL::Origin& /*o*/) : cx(0), cy(0), cz(0) { }
  My_3D_Point(const FT x, const FT y, const FT z) : cx(x), cy(y), cz(z) { }

  FT x() const { return cx; }
  FT y() const { return cy; }
  FT z() const { return cz; }

private:
  FT cx, cy, cz;
};

template <typename K, typename VPM>
struct Locate_with_AABB_tree_Tester<K, VPM, 3> // 3D
{
  template <typename G>
  void test(const G& g, const VPM vpm, CGAL::Random& rnd) const
  {
    std::cout << "  test locate_with_AABB_tree (3D)..." << std::endl;

    typedef typename boost::property_traits<VPM>::reference                    Point_reference;

    typedef typename K::FT                                                     FT;
    typedef typename K::Ray_3                                                  Ray_3;

    typedef typename boost::graph_traits<G>::vertex_descriptor                 vertex_descriptor;
    typedef typename boost::graph_traits<G>::halfedge_descriptor               halfedge_descriptor;
    typedef typename boost::graph_traits<G>::face_descriptor                   face_descriptor;

    typedef typename PMP::Face_location<G, FT>                                 Face_location;

    face_descriptor f = CGAL::internal::random_face_in_mesh(g, rnd);
    halfedge_descriptor h = halfedge(f, g);
    vertex_descriptor v = target(h, g);

    // ---------------------------------------------------------------------------
    typedef CGAL::AABB_face_graph_triangle_primitive<G, VPM>                   AABB_face_graph_primitive;
    typedef CGAL::AABB_traits<K, AABB_face_graph_primitive>                    AABB_face_graph_traits;

    typedef typename K::Point_3                            Point_3;
    static_assert(std::is_same<typename AABB_face_graph_traits::Point_3, Point_3>::value);

    CGAL::AABB_tree<AABB_face_graph_traits> tree_a;
    Point_reference p3_a = get(vpm, v);

    // below tests the case where the value type of the VPM is not Kernel::Point_3
    typedef My_3D_Point<K>                                                     Custom_point;
    typedef std::map<vertex_descriptor, Custom_point>                          Custom_map;
    typedef boost::associative_property_map<Custom_map>                        Custom_VPM;
    typedef PMP::internal::Point_to_Point_3_VPM<G, Custom_VPM>                 WrappedVPM;
    typedef CGAL::AABB_face_graph_triangle_primitive<G, WrappedVPM>            AABB_face_graph_primitive_with_WVPM;
    typedef CGAL::AABB_traits<K, AABB_face_graph_primitive_with_WVPM>          AABB_face_graph_traits_with_WVPM;

    CGAL::AABB_tree<AABB_face_graph_traits_with_WVPM> tree_b;
    Custom_map custom_map;
    for(vertex_descriptor vd : vertices(g))
    {
      const Point_reference p = get(vpm, vd);
      custom_map[vd] = Custom_point(p.x(), p.y(), p.z());
    }

    Custom_VPM custom_vpm(custom_map);
    WrappedVPM custom_vpm_3D(custom_vpm);
    // ---------------------------------------------------------------------------

    PMP::build_AABB_tree(g, tree_a); // just for the API
    assert(tree_a.size() == num_faces(g));

    PMP::build_AABB_tree(g, tree_a, CGAL::parameters::vertex_point_map(vpm));
    PMP::build_AABB_tree(g, tree_b, CGAL::parameters::vertex_point_map(custom_vpm_3D));
    assert(tree_b.size() == num_faces(g));

    Face_location loc = PMP::locate_with_AABB_tree(p3_a, tree_a, g, CGAL::parameters::vertex_point_map(vpm));
    if (std::is_same<K, EPECK>()) {
      assert(is_equal(loc.second[CGAL::vertex_index_in_face(v, loc.first, g)], FT(1)));
      assert(is_equal(loc.second[(CGAL::vertex_index_in_face(v, loc.first, g) + 1) % 3], FT(0)));
      assert(is_equal(loc.second[(CGAL::vertex_index_in_face(v, loc.first, g) + 2) % 3], FT(0)));
      assert(is_equal(CGAL::squared_distance(PMP::construct_point(loc, g), p3_a), FT(0)));
    }

    loc = PMP::locate_with_AABB_tree(p3_a, tree_a, g, CGAL::parameters::vertex_point_map(vpm));
    if (std::is_same<K, EPECK>()) {
      assert(is_equal(CGAL::squared_distance(PMP::construct_point(loc, g), p3_a), FT(0)));
    }

    // ---------------------------------------------------------------------------
    loc = PMP::locate(p3_a, g, CGAL::parameters::snapping_tolerance(1e-7));
    if (std::is_same<K, EPECK>()) {
      assert(is_equal(CGAL::squared_distance(PMP::construct_point(loc, g), p3_a), FT(0)));
      assert(PMP::is_in_face(loc, g));
    }

    loc = PMP::locate_with_AABB_tree(CGAL::ORIGIN, tree_b, g, CGAL::parameters::vertex_point_map(custom_vpm_3D));
    if (std::is_same<K, EPECK>()) {
      assert(PMP::is_in_face(loc, g));
    }

    // Doesn't necessarily have to wrap with a P_to_P3: it can be done automatically internally
    loc = PMP::locate(CGAL::ORIGIN, g, CGAL::parameters::vertex_point_map(custom_vpm));
    if (std::is_same<K, EPECK>()) {
      assert(PMP::is_in_face(loc, g));
    }

    // ---------------------------------------------------------------------------
    Ray_3 r3 = random_3D_ray<CGAL::AABB_tree<AABB_face_graph_traits_with_WVPM> >(tree_b, rnd);
    loc = PMP::locate_with_AABB_tree(r3, tree_b, g, CGAL::parameters::vertex_point_map(custom_vpm_3D));
  }
};


template<typename K, typename G, typename VPM>
void test_locate(const G& g,
                 const VPM vpm,
                 CGAL::Random& rnd)
{
  assert(num_vertices(g) != 0 && num_faces(g) != 0);

  test_snappers<K>(g);
  test_constructions<K>(g, vpm, rnd);
  test_random_entities<K>(g, rnd);
  test_helpers<K>(g, rnd);
  test_predicates<K>(g, rnd);
  test_locate_in_face<K>(g, vpm, rnd);

  // This test has slight syntax changes between 2D and 3D (e.g. testing ray_2 in 3D makes no sense)
  Locate_with_AABB_tree_Tester<K, VPM> AABB_tester;
  AABB_tester.test(g, vpm, rnd);
}

template<typename K, typename G>
void test_locate(const G& g, CGAL::Random& rnd)
{
  return test_locate<K>(g, CGAL::get_const_property_map(boost::vertex_point, g), rnd);
}

template<typename K>
void test_2D_triangulation(const std::string fname, CGAL::Random& rnd)
{
  typedef CGAL::Regular_triangulation_2<K>                    RT;

  std::cout << "Testing Regular_triangulation_2 " << fname;

  std::ifstream in(fname);
  assert(in.good());

  RT tr;
  double x, y;
  while(in >> x >> y)
    tr.insert(typename RT::Point(x, y));

  std::ofstream out("triangulation.off");
  out << "OFF\n";
  out << tr.number_of_vertices() << " " << std::distance(tr.finite_faces_begin(), tr.finite_faces_end()) << " 0\n";

  std::size_t counter = 0;
  std::map<typename RT::Point_2, std::size_t> ids;
  for(const auto& v : CGAL::make_range(tr.finite_vertices_begin(), tr.finite_vertices_end()))
  {
    out << v.point().point() << " 0\n";
    if(ids.insert(std::make_pair(v.point().point(), counter)).second)
      ++counter;
  }

  for(const auto& fd : CGAL::make_range(tr.finite_faces_begin(), tr.finite_faces_end()))
  {
    out << "3 " << ids[fd.vertex(0)->point().point()] << " " << ids[fd.vertex(1)->point().point()] << " " << ids[fd.vertex(2)->point().point()] << "\n";
  }

  out.close();

  std::cout << " (" << tr.number_of_vertices() << " vertices)..." << std::endl;
  std::cout << "Kernel: " << typeid(K()).name() << std::endl;

  test_locate<K>(tr, rnd);
}

template<typename K>
void test_2D_surface_mesh(const std::string fname, CGAL::Random& rnd)
{
  typedef typename K::Point_2                                 Point;
  typedef CGAL::Surface_mesh<Point>                           Mesh;

  std::cout << "Testing Surface_mesh " << fname << "..." << std::endl;
  std::cout << "Kernel: " << typeid(K()).name() << std::endl;

  std::ifstream input(fname);
  assert(input.good());

  Mesh tm;
  if(!input || !(input >> tm))
  {
    std::cerr << "Error: cannot read file.";
    return;
  }

  test_locate<K>(tm, rnd);
}

template<typename K>
void test_surface_mesh_3D(const std::string fname, CGAL::Random& rnd)
{
  typedef typename K::Point_3                                 Point;
  typedef CGAL::Surface_mesh<Point>                           Mesh;

  std::cout << "Testing (3D) Surface_mesh " << fname << "..." << std::endl;
  std::cout << "Kernel: " << typeid(K()).name() << std::endl;

  std::ifstream input(fname);
  Mesh tm;
  if(!input || !(input >> tm))
  {
    std::cerr << "Error: cannot read file.";
    return;
  }

  typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::const_type  VertexPointMap;
  VertexPointMap vpm = CGAL::get_const_property_map(boost::vertex_point, tm);

  test_locate<K>(tm, vpm, rnd);
}

template<typename K>
void test_surface_mesh_projection(const std::string fname, CGAL::Random& rnd)
{
  typedef typename K::Point_3                                       Point;
  typedef CGAL::Surface_mesh<Point>                                 Mesh;
  typedef typename K::Point_2                                       Projected_point;

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor     vertex_descriptor;

  std::cout << "Testing Projected Surface_mesh " << fname << "..." << std::endl;
  std::cout << "Kernel: " << typeid(K()).name() << std::endl;

  std::ifstream input(fname);
  Mesh tm;
  if(!input || !(input >> tm))
  {
    std::cerr << "Error: cannot read file.";
    return;
  }

  const auto& proj_vpm = tm.template add_property_map<typename Mesh::Vertex_index,
                                                      Projected_point>("P2", Projected_point()).first;

  for(vertex_descriptor v : vertices(tm))
  {
    const Point& p = tm.point(v);
    put(proj_vpm, v, Projected_point(p.x(), p.y()));
  }

  test_locate<K>(tm, proj_vpm, rnd);
}

template<typename K>
void test_polyhedron(const std::string fname, CGAL::Random& rnd)
{
  typedef CGAL::Polyhedron_3<K>                               Polyhedron;

  std::cout << "Testing Polyhedron_3 " << fname << "..." << std::endl;
  std::cout << "Kernel: " << typeid(K()).name() << std::endl;

  std::ifstream input(fname);
  Polyhedron poly;
  if(!input || !(input >> poly))
  {
    std::cerr << "Error: cannot read file.";
    return;
  }

  test_locate<K>(poly, rnd);
}

template <typename K>
void test(CGAL::Random& rnd)
{
  test_2D_triangulation<K>("data/stair.xy", rnd);
//  test_2D_surface_mesh<K>("data/blobby_2D.off", rnd); // temporarily disabled, until Surface_mesh's IO is "fixed"
  test_surface_mesh_3D<K>("meshes/mech-holes-shark.off", rnd);
  test_surface_mesh_projection<K>("data/unit-grid.off", rnd);
  test_polyhedron<K>("data-coref/elephant_split_2.off", rnd);
}

int main()
{
  CGAL::Set_ieee_double_precision pfr;

  std::cout.precision(17);
  std::cout << std::fixed;

//  CGAL::Random rnd(1557332474); // if needed to debug with a fixed seed
  CGAL::Random rnd(CGAL::get_default_random());

  std::cout << "The seed is " << rnd.get_seed() << std::endl;

  test<EPICK>(rnd);
  test<Exact_kernel>(rnd);

  return EXIT_SUCCESS;
}
