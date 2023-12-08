#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/use.h>
#include "test_Prefix.h"

#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <map>
#include <memory>
#include <utility>
#include <cassert>

typedef std::unordered_set<std::size_t> id_map;

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename Graph>
void test_halfedge_around_vertex_iterator(const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::face_descriptor g_face_descriptor;
  typedef CGAL::Face_filtered_graph<Graph> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::unordered_map<g_face_descriptor, std::size_t> map(num_faces(g));
  PMP::connected_components(g, boost::make_assoc_property_map(map));

  Adapter fg(g, 0, boost::make_assoc_property_map(map));
  typename boost::graph_traits<Adapter >::vertex_iterator vit, vend;
  for(boost::tie(vit, vend) = vertices(fg); vit != vend; ++vit) {
    halfedge_around_target_iterator havit, havend;
    for(boost::tie(havit, havend) = CGAL::halfedges_around_target(halfedge(*vit, fg), fg);
        havit != havend; ++havit) {
      assert(target(*havit, fg) == *vit);

      // check if we are really moving clockwise
      halfedge_around_target_iterator step = std::next(havit);
      if(step != havend) {
        halfedge_descriptor stepd = *step;
        assert(stepd == opposite(next(*havit, fg), fg));
      }
    }
  }
}

template <typename Graph>
void test_halfedge_around_face_iterator(const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::face_descriptor g_face_descriptor;
  typedef CGAL::Face_filtered_graph<Graph> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  PMP::connected_components(g, boost::make_assoc_property_map(map));
  Adapter fg(g, 0, boost::make_assoc_property_map(map));

  face_iterator fit, fend;
  for(boost::tie(fit, fend) = faces(fg); fit != fend; ++fit) {
    halfedge_around_face_iterator hafit, hafend;
    boost::tie(hafit, hafend) = CGAL::halfedges_around_face(halfedge(*fit, fg), fg);
    assert(std::distance(hafit, hafend) != 0);
    for(boost::tie(hafit, hafend) = CGAL::halfedges_around_face(halfedge(*fit, fg), fg); hafit != hafend; ++hafit) {
      assert(face(*hafit, fg) == *fit);
    }
  }
}

template<typename Graph>
void test_edge_iterators(const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::face_descriptor g_face_descriptor;
  typedef CGAL::Face_filtered_graph<Graph> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  PMP::connected_components(g, boost::make_assoc_property_map(map));
  Adapter fg(g, 0, boost::make_assoc_property_map(map));

  // do we iterate as many as that?
  edge_iterator eb, ee;
  boost::tie(eb, ee) = edges(fg);
  assert(static_cast<edges_size_type>(std::distance(eb, ee)) == num_edges(g));
  id_map ids;
  unsigned int count = 0;
  for(boost::tie(eb, ee) = edges(fg); eb != ee; ++eb) {
    edge_descriptor e = *eb;
    std::pair<id_map::iterator, bool> r = ids.insert(get(boost::edge_index, g, e));
    // unique?
    assert(r.second);
    ++count;
  }
  assert(count == num_edges(fg));
}

template<typename Graph>
void test_vertex_iterators(Graph& g)
{
  typedef typename boost::graph_traits<Graph>::face_descriptor g_face_descriptor;
  typedef CGAL::Face_filtered_graph<Graph> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  PMP::connected_components(g, boost::make_assoc_property_map(map));
  Adapter fg(g, 0, boost::make_assoc_property_map(map));
  vertex_iterator vb, ve;
  std::size_t count = 0;
  for(boost::tie(vb, ve) = vertices(fg); vb != ve; ++vb){
    ++count;
  }

  assert(count == num_vertices(fg));

  // check that the iterators reach uniques
  id_map ids;

  count = 0;
  for(boost::tie(vb, ve) = vertices(fg); vb != ve; ++vb) {
    std::pair<id_map::iterator, bool> r = ids.insert(get(boost::vertex_index, g, *vb));
    assert(r.second);
    ++count;
  }
  assert(count == num_vertices(fg));
}


template<typename Graph>
void test_out_edges(const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::face_descriptor g_face_descriptor;
  typedef CGAL::Face_filtered_graph<Graph> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  PMP::connected_components(g, boost::make_assoc_property_map(map));
  Adapter fg(g, 0, boost::make_assoc_property_map(map));

  vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = vertices(fg); vb != ve; ++vb) {
    id_map v_ids;

    vertex_descriptor around = *vb;
    out_edge_iterator oeb, oee;
    for(boost::tie(oeb, oee) = out_edges(*vb, fg); oeb != oee; ++oeb) {
      vertex_descriptor t = target(*oeb, fg);
      vertex_descriptor s = source(*oeb, fg);
      assert(s != t);
      assert(s == around);
      assert(t != around);
      std::pair<id_map::iterator, bool> r =
          v_ids.insert(get(boost::vertex_index, g, target(*oeb, fg)));
      assert(r.second);
    }
  }
}

template<typename Graph>
void test_in_edges(const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::face_descriptor g_face_descriptor;
  typedef CGAL::Face_filtered_graph<Graph> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  PMP::connected_components(g, boost::make_assoc_property_map(map));
  Adapter fg(g, 0, boost::make_assoc_property_map(map));

  vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = vertices(fg); vb != ve; ++vb) {
    id_map v_ids;
    vertex_descriptor around = *vb;
    in_edge_iterator ieb, iee;
    for(boost::tie(ieb, iee) = in_edges(*vb, fg); ieb != iee; ++ieb) {
      vertex_descriptor t = target(*ieb, fg);
      vertex_descriptor s = source(*ieb, fg);
      assert(t == around);
      assert(s != around);
      std::pair<id_map::iterator, bool> r =
          v_ids.insert(get(boost::vertex_index, g, source(*ieb, fg)));
      assert(r.second);
    }
  }
}

template<typename Graph>
void test_in_out_edges(const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::face_descriptor g_face_descriptor;
  typedef CGAL::Face_filtered_graph<Graph> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  PMP::connected_components(g, boost::make_assoc_property_map(map));
  Adapter fg(g, 0, boost::make_assoc_property_map(map));

  // check that the sets of in out edges are the same
  vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = vertices(fg); vb != ve; ++vb) {
    id_map v_ids;
    std::vector<vertex_descriptor> in, out;
    in_edge_iterator ieb, iee;
    for(boost::tie(ieb, iee) = in_edges(*vb, fg); ieb != iee; ++ieb) {
      std::pair<id_map::iterator, bool> r =
          v_ids.insert(get(boost::vertex_index, g, source(*ieb, fg)));
      assert(r.second);
      in.push_back(source(*ieb, fg));
    }
    out_edge_iterator oeb, oee;
    for(boost::tie(oeb, oee) = out_edges(*vb, fg); oeb != oee; ++oeb) {
      std::pair<id_map::iterator, bool> r =
          v_ids.insert(get(boost::vertex_index, g, target(*oeb, fg)));
      // insertion must fail
      assert(!r.second);
      out.push_back(target(*oeb, fg));
    }
    // did we walk the vertices in the same order?
    assert(in.size() == out.size());
    assert(std::equal(in.begin(), in.end(), out.begin()));
    assert(in.size() == in_degree(*vb, fg));
    assert(out.size() == out_degree(*vb, fg));
    assert(in.size() == degree(*vb, fg));
    assert(degree(*vb, fg) == in_degree(*vb, fg));
    assert(degree(*vb, fg) == out_degree(*vb, fg));
  }
}

// check that every edge can be found through edge(u, v, g)
template<typename Graph>
void test_edge_find(const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::face_descriptor g_face_descriptor;
  typedef CGAL::Face_filtered_graph<Graph> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  PMP::connected_components(g, boost::make_assoc_property_map(map));
  Adapter fg(g, 0, boost::make_assoc_property_map(map));
  typedef std::pair<edge_descriptor, bool>   ret;

  edge_iterator eb, ee;
  for(boost::tie(eb, ee) = edges(fg); eb != ee; ++eb) {
    vertex_descriptor s = source(*eb, fg);
    vertex_descriptor t = target(*eb, fg);
    ret found = edge(s, t, fg);
    ret found2 = edge(t, s, fg);
    assert(found.second);
    assert(found2.second);
    assert(found.first == *eb);
    assert(found2.first == *eb);
  }
}

template<typename Graph>
void test_faces(const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::face_descriptor g_face_descriptor;
  typedef CGAL::Face_filtered_graph<Graph> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  PMP::connected_components(g, boost::make_assoc_property_map(map));
  Adapter fg(g, 0, boost::make_assoc_property_map(map));

  unsigned int count = 0;
  face_iterator fb, fe;
  for(boost::tie(fb, fe) = faces(fg); fb != fe; ++fb) {
    ++count;
    // reverse look-up
    halfedge_descriptor assoc = halfedge(*fb, fg);
    assert(face(assoc, fg) == *fb);
    // check the enclosure
    halfedge_around_face_iterator encb, ence;
    for(boost::tie(encb, ence) = CGAL::halfedges_around_face(halfedge(*fb, fg), fg); encb != ence; ++encb) {
      assert(face(*encb, fg) == *fb);
    }
  }
  assert(count == num_faces(fg));
}

template<typename Graph>
void test_index_property_maps(const Graph& g)
{
  typedef CGAL::Face_filtered_graph<Graph> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);

  typedef typename boost::graph_traits<Graph>::face_descriptor g_face_descriptor;
  std::map<g_face_descriptor, std::size_t> map;
  PMP::connected_components(g, boost::make_assoc_property_map(map), CGAL::parameters::default_values());
  Adapter fg(g, -1, boost::make_assoc_property_map(map));
  assert(is_empty(fg));

  // Non const
  typedef typename CGAL::GetInitializedVertexIndexMap<Adapter>::type VIMap;
  VIMap vim = CGAL::get_initialized_vertex_index_map(fg);
  assert(CGAL::BGL::internal::is_index_map_valid(vim, num_vertices(fg), vertices(fg)));

  typedef typename CGAL::GetInitializedHalfedgeIndexMap<Adapter>::type HIMap;
  HIMap him = CGAL::get_initialized_halfedge_index_map(fg);
  assert(CGAL::BGL::internal::is_index_map_valid(him, num_halfedges(fg), halfedges(fg)));

  typedef typename CGAL::GetInitializedFaceIndexMap<Adapter>::type FIMap;
  FIMap fim = CGAL::get_initialized_face_index_map(fg);
  assert(CGAL::BGL::internal::is_index_map_valid(fim, num_faces(fg), faces(fg)));

  fg.set_selected_faces(0, boost::make_assoc_property_map(map));
  assert(!is_empty(fg));

  assert(CGAL::BGL::internal::is_index_map_valid(vim, num_vertices(fg), vertices(fg)));
  assert(CGAL::BGL::internal::is_index_map_valid(him, num_halfedges(fg), halfedges(fg)));
  assert(CGAL::BGL::internal::is_index_map_valid(fim, num_faces(fg), faces(fg)));

  // Const
  const Adapter cfg(g, 0, boost::make_assoc_property_map(map));
  typedef typename CGAL::GetInitializedVertexIndexMap<Adapter>::const_type CVIMap;
  CVIMap cvim = CGAL::get_initialized_vertex_index_map(cfg);
  assert(CGAL::BGL::internal::is_index_map_valid(cvim, num_vertices(cfg), vertices(cfg)));

  typedef typename CGAL::GetInitializedHalfedgeIndexMap<Adapter>::const_type CHIMap;
  CHIMap chim = CGAL::get_initialized_halfedge_index_map(cfg);
  assert(CGAL::BGL::internal::is_index_map_valid(chim, num_halfedges(cfg), halfedges(cfg)));

  typedef typename CGAL::GetInitializedFaceIndexMap<Adapter>::const_type CFIMap;
  CFIMap cfim = CGAL::get_initialized_face_index_map(cfg);
  assert(CGAL::BGL::internal::is_index_map_valid(cfim, num_faces(cfg), faces(cfg)));
}

template<typename Graph>
void test_constructors(const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::face_descriptor g_face_descriptor;
  typedef CGAL::Face_filtered_graph<Graph> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);

  Adapter fg0(g);
  Adapter fg1(g, CGAL::parameters::default_values());

  std::map<g_face_descriptor, std::size_t> map;
  PMP::connected_components(g, boost::make_assoc_property_map(map));
  Adapter fg(g, 0, boost::make_assoc_property_map(map));
  assert(fg.is_selection_valid());
  assert(CGAL::is_valid_polygon_mesh(fg));

  // patch ID type that is not a fundamental type
  std::map<g_face_descriptor, CGAL::IO::Color> colors;
  for(auto f : faces(fg))
    colors[f] = CGAL::IO::blue();
  auto color_map = boost::make_assoc_property_map(colors);
  Adapter color_fg(g, CGAL::IO::blue(), color_map);
  assert(color_fg.is_selection_valid());
  assert(CGAL::is_valid_polygon_mesh(color_fg));
  assert(num_faces(g) == num_faces(color_fg));

  Adapter color_fg2(g, CGAL::IO::red(), color_map, CGAL::parameters::default_values());
  assert(color_fg2.is_selection_valid());
  assert(CGAL::is_valid_polygon_mesh(color_fg2));
  assert(num_faces(color_fg2) == 0);

  // range of non-fundamental ID types
  std::vector<CGAL::IO::Color> selected_colors { CGAL::IO::blue(), CGAL::IO::red() };
  Adapter color_fg3(g, selected_colors, color_map);
  assert(color_fg3.is_selection_valid());
  assert(CGAL::is_valid_polygon_mesh(color_fg3));
  assert(num_faces(g) == num_faces(color_fg3));
}

template <typename Graph>
void
test_graph_range(const std::vector<Graph>& graphs)
{
  for(Graph p : graphs)
  {
    test_constructors(p);
    test_vertex_iterators(p);
    test_out_edges(p);
    test_in_edges(p);
    test_in_out_edges(p);
    test_edge_find(p);
    test_faces(p);
    test_edge_iterators(p);
    test_halfedge_around_face_iterator(p);
    test_halfedge_around_vertex_iterator(p);
    test_index_property_maps(p);
  }
}


typedef SM::Point Point_3;

template<class Mesh, typename VertexPointPMap>
struct Constraint
{
  typedef typename boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;
  typedef boost::readable_property_map_tag      category;
  typedef bool                                  value_type;
  typedef bool                                  reference;
  typedef edge_descriptor                       key_type;

  Constraint()
    :g_(NULL)
  {}

  Constraint(Mesh& g, VertexPointPMap vpp)
    : g_(&g), vppmap(vpp)
  {}

  bool operator[](edge_descriptor e) const
  {
    const Mesh& g = *g_;
    if(
       ((boost::get(vppmap, target(e, g)) == Point_3(1,1,1) ||
         boost::get(vppmap, source(e, g)) == Point_3(1,1,1)) &&
        (boost::get(vppmap, target(e, g)) == Point_3(0,0,0) ||
         boost::get(vppmap, source(e, g)) == Point_3(0,0,0)))
       ||
       ((boost::get(vppmap, target(e, g)) == Point_3(1,1,1) ||
         boost::get(vppmap, source(e, g)) == Point_3(1,1,1)) &&
        (boost::get(vppmap, target(e, g)) == Point_3(0,0,1) ||
         boost::get(vppmap, source(e, g)) == Point_3(0,0,1)))
       ||
       ((boost::get(vppmap, target(e, g)) == Point_3(0,0,1) ||
         boost::get(vppmap, source(e, g)) == Point_3(0,0,1)) &&
        (boost::get(vppmap, target(e, g)) == Point_3(0,0,0) ||
         boost::get(vppmap, source(e, g)) == Point_3(0,0,0)))
       ||
       ((boost::get(vppmap, target(e, g)) == Point_3(1,0,1) ||
         boost::get(vppmap, source(e, g)) == Point_3(1,0,1)) &&
        (boost::get(vppmap, target(e, g)) == Point_3(0,0,0) ||
         boost::get(vppmap, source(e, g)) == Point_3(0,0,0)))
       ||
       ((boost::get(vppmap, target(e, g)) == Point_3(1,1,1) ||
         boost::get(vppmap, source(e, g)) == Point_3(1,1,1)) &&
        (boost::get(vppmap, target(e, g)) == Point_3(1,0,1) ||
         boost::get(vppmap, source(e, g)) == Point_3(1,0,1)))
       )
      return true;
    else
      return false;
  }

  friend inline value_type get(const Constraint& m, const key_type k) { return m[k]; }

  const Mesh* g_;
  VertexPointPMap vppmap;
};

template<class Mesh, class FCCMAP, class Adapter>
void test_mesh(Adapter fga)
{
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  //check that there is the right number of simplices in fga
  assert(CGAL::is_valid_polygon_mesh(fga));
  assert(num_faces(fga) == 2);
  assert(num_edges(fga) == 5);
  assert(num_halfedges(fga) == 10);
  assert(num_vertices(fga) == 4);
  halfedge_descriptor h = halfedge(*faces(fga).first, fga);
  vertex_descriptor v = source(h, fga);
  //check that next() works inside the patch
  assert(next(next(next(h, fga), fga), fga) == h);
  //check that next() works on bordure of the patch
  h = opposite(h, fga);
  assert(next(next(next(next(h, fga), fga), fga), fga) == h);
  //check that prev() works inside the patch
  h = halfedge(*faces(fga).first, fga);
  assert(prev(prev(prev(h, fga), fga), fga) == h);
  //check that prev() works on bordure of the patch
  h = opposite(h, fga);
  assert(prev(prev(prev(prev(h, fga), fga), fga), fga) == h);
  //check degree
  assert(degree(v, fga) == 3);
  //check in_edges and out_edges
  assert(std::distance(in_edges(v, fga).first ,in_edges(v, fga).second) == 3);
  assert(std::distance(out_edges(v, fga).first ,out_edges(v, fga).second) == 3);

  Mesh copy;
  CGAL::copy_face_graph(fga, copy);
}

template <typename PolygonMesh>
void merge_vertices(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v_keep,
                    typename boost::graph_traits<PolygonMesh>::vertex_descriptor v_rm,
                    std::vector<typename boost::graph_traits<PolygonMesh>::face_descriptor>& incident_faces,
                    PolygonMesh& mesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor oh = halfedge(v_keep, mesh), done = oh;
  do
  {
    incident_faces.push_back(face(oh, mesh));
    oh = opposite(next(oh, mesh), mesh);
  }
  while(oh != done);

  halfedge_descriptor h = halfedge(v_rm, mesh);
  halfedge_descriptor start = h;
  do
  {
    set_target(h, v_keep, mesh);
    incident_faces.push_back(face(h, mesh));
    h = opposite(next(h, mesh), mesh);
  }
  while(h != start);

  remove_vertex(v_rm, mesh);
}

void test_invalid_selections()
{
  // this creates a non-manifold (pinched) vertex
  SM mesh;
  read_a_mesh(mesh, "data/7_faces_triangle.off");

  std::vector<SM::Face_index> face_range;
  face_range.push_back(SM::Face_index(1));
  face_range.push_back(SM::Face_index(2));
  face_range.push_back(SM::Face_index(3));

  CGAL::Face_filtered_graph<SM> pinched_fg(mesh, face_range);
  assert(pinched_fg.is_selection_valid());

  // this creates a non-manifold vertex (multiple umbrellas)
  clear(mesh);
  read_a_mesh(mesh, "data/genus3.off");
  assert(is_valid_polygon_mesh(mesh));

  face_range.clear();
  merge_vertices(SM::Vertex_index(1337), SM::Vertex_index(87), face_range, mesh);

  CGAL::Face_filtered_graph<SM> many_umbrellas_fg(mesh, face_range);
  assert(!many_umbrellas_fg.is_selection_valid());
}

void test_SM_tetrahedron()
{
  // Make a tetrahedron and test the adapter for a patch that only contains 2 faces
  typedef CGAL::Face_filtered_graph<SM> SM_Adapter;
  typedef SM::Property_map<boost::graph_traits<SM>::face_descriptor , std::size_t> SM_FCCMap;
  auto sm = std::make_unique<SM>();
  CGAL::make_tetrahedron(Point_3(1,1,1), Point_3(0,0,0), Point_3(0,0,1), Point_3(1,0,1), *sm);

  SM_FCCMap fccmap = sm->add_property_map<boost::graph_traits<SM>::face_descriptor, std::size_t>("f:CC").first;
  SM::Property_map<boost::graph_traits<SM>::vertex_descriptor, SM::Point> positions = sm->points();
  CGAL::Polygon_mesh_processing::connected_components(
    *sm, fccmap, CGAL::parameters::edge_is_constrained_map(Constraint<SM, SM::Property_map<boost::graph_traits<SM>::vertex_descriptor,
                                                           SM::Point> >(*sm, positions)));

  std::unordered_set<long unsigned int> pids;
  pids.insert(0);
  pids.insert(2);
  SM_Adapter sm_adapter(*sm, pids, fccmap);
  test_mesh<SM,SM_FCCMap, SM_Adapter>(sm_adapter);
}

void test_Polyhedron_tetrahedron()
{
  typedef boost::graph_traits<Polyhedron> PolyTraits;
  typedef boost::property_map<Polyhedron, boost::vertex_point_t>::const_type VPMap;
  typedef PolyTraits::face_descriptor poly_face_descriptor;
  typedef boost::associative_property_map<std::map<poly_face_descriptor, PolyTraits::faces_size_type> > FCMap;
  typedef boost::property_map<Polyhedron, CGAL::face_external_index_t>::const_type FIMap;
  typedef boost::property_map<Polyhedron, CGAL::vertex_external_index_t>::const_type VIMap;
  typedef boost::property_map<Polyhedron, CGAL::halfedge_external_index_t>::const_type HIMap;
  typedef CGAL::Face_filtered_graph<Polyhedron, FIMap, VIMap, HIMap> Poly_Adapter;

  auto poly = std::make_unique<Polyhedron>();
  CGAL::make_tetrahedron(Point_3(1,1,1), Point_3(0,0,0), Point_3(0,0,1), Point_3(1,0,1), *poly);

  FIMap poly_fimap = get(CGAL::face_external_index, *poly);
  VIMap poly_vimap = get(CGAL::vertex_external_index, *poly);
  HIMap poly_himap = get(CGAL::halfedge_external_index, *poly);
  std::map<poly_face_descriptor, PolyTraits::faces_size_type> fc_map;
  FCMap poly_fccmap(fc_map);

  std::unordered_set<long unsigned int> pids;
  pids.insert(0);
  pids.insert(2);
  VPMap vpmap = get(boost::vertex_point, *poly);
  PMP::connected_components(*poly, poly_fccmap,
                            CGAL::parameters::edge_is_constrained_map(Constraint<Polyhedron, VPMap >(*poly, vpmap))
                                             .face_index_map(poly_fimap));
  Poly_Adapter poly_adapter(*poly, pids, poly_fccmap,
                            CGAL::parameters::face_index_map(poly_fimap)
                                             .vertex_index_map(poly_vimap)
                                             .halfedge_index_map(poly_himap));
  test_mesh<Polyhedron, FCMap, Poly_Adapter>(poly_adapter);
}

int main(int, char**)
{
  test_graph_range(poly_data());
#ifdef CGAL_USE_SURFACE_MESH
  test_graph_range(sm_data());
#endif
#ifdef CGAL_USE_OPENMESH
  test_graph_range(omesh_data());
#endif

  test_invalid_selections();

  test_SM_tetrahedron();
  test_Polyhedron_tetrahedron();

  std::cout << "Done" << std::endl;
  return EXIT_SUCCESS;
}
