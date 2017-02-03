#include <CGAL/boost/graph/Connected_components_graph.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include "test_Prefix.h"
#include <boost/numeric/conversion/cast.hpp>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <CGAL/use.h>

typedef boost::unordered_set<std::size_t> id_map;

template <typename Graph>
void test_halfedge_around_vertex_iterator(const  Graph& g)
{
  typedef typename boost::graph_traits<Graph>::face_descriptor g_face_descriptor;
  typedef boost::associative_property_map< boost::unordered_map< g_face_descriptor, std::size_t > >FCMap;
  typedef CGAL::Connected_components_graph<Graph, FCMap> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  boost::unordered_map<g_face_descriptor, std::size_t> map(CGAL::num_faces(g));
  CGAL::Polygon_mesh_processing::connected_components(g, boost::make_assoc_property_map(map), CGAL::Polygon_mesh_processing::parameters::all_default());


  Adapter fg(g, boost::make_assoc_property_map(map), 0);
  typename boost::graph_traits<Adapter >::vertex_iterator vit, vend;
  for(boost::tie(vit, vend) = vertices(fg); vit != vend; ++vit) {
    halfedge_around_target_iterator havit, havend;
    for(boost::tie(havit, havend) = CGAL::halfedges_around_target(halfedge(*vit, fg), fg);
        havit != havend; ++havit) {
      assert(target(*havit, fg) == *vit);

      // check if we are really moving clockwise
      halfedge_around_target_iterator step = boost::next(havit);
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
  typedef boost::associative_property_map< std::map< g_face_descriptor, std::size_t > >FCMap;
  typedef CGAL::Connected_components_graph<Graph, FCMap> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  CGAL::Polygon_mesh_processing::connected_components(g, boost::make_assoc_property_map(map), CGAL::Polygon_mesh_processing::parameters::all_default());
  Adapter fg(g, boost::make_assoc_property_map(map), 0);
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
  typedef boost::associative_property_map< std::map< g_face_descriptor, std::size_t > >FCMap;
  typedef CGAL::Connected_components_graph<Graph, FCMap> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  CGAL::Polygon_mesh_processing::connected_components(g, boost::make_assoc_property_map(map), CGAL::Polygon_mesh_processing::parameters::all_default());
  Adapter fg(g, boost::make_assoc_property_map(map), 0);

  // do we iterate as many as that?
  edge_iterator eb, ee;
  boost::tie(eb, ee) = edges(fg);
  assert(boost::numeric_cast<edges_size_type>(std::distance(eb, ee)) == num_edges(g));
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
  typedef boost::associative_property_map< std::map< g_face_descriptor, std::size_t > >FCMap;
  typedef CGAL::Connected_components_graph<Graph, FCMap> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  CGAL::Polygon_mesh_processing::connected_components(g, boost::make_assoc_property_map(map), CGAL::Polygon_mesh_processing::parameters::all_default());
  Adapter fg(g, boost::make_assoc_property_map(map), 0);
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
  typedef boost::associative_property_map< std::map< g_face_descriptor, std::size_t > >FCMap;
  typedef CGAL::Connected_components_graph<Graph, FCMap> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  CGAL::Polygon_mesh_processing::connected_components(g, boost::make_assoc_property_map(map), CGAL::Polygon_mesh_processing::parameters::all_default());
  Adapter fg(g, boost::make_assoc_property_map(map), 0);

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
  typedef boost::associative_property_map< std::map< g_face_descriptor, std::size_t > >FCMap;
  typedef CGAL::Connected_components_graph<Graph, FCMap> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  CGAL::Polygon_mesh_processing::connected_components(g, boost::make_assoc_property_map(map), CGAL::Polygon_mesh_processing::parameters::all_default());
  Adapter fg(g, boost::make_assoc_property_map(map), 0);

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
  typedef boost::associative_property_map< std::map< g_face_descriptor, std::size_t > >FCMap;
  typedef CGAL::Connected_components_graph<Graph, FCMap> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  CGAL::Polygon_mesh_processing::connected_components(g, boost::make_assoc_property_map(map), CGAL::Polygon_mesh_processing::parameters::all_default());
  Adapter fg(g, boost::make_assoc_property_map(map), 0);

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
  typedef boost::associative_property_map< std::map< g_face_descriptor, std::size_t > >FCMap;
  typedef CGAL::Connected_components_graph<Graph, FCMap> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  CGAL::Polygon_mesh_processing::connected_components(g, boost::make_assoc_property_map(map), CGAL::Polygon_mesh_processing::parameters::all_default());
  Adapter fg(g, boost::make_assoc_property_map(map), 0);
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
  typedef boost::associative_property_map< std::map< g_face_descriptor, std::size_t > >FCMap;
  typedef CGAL::Connected_components_graph<Graph, FCMap> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  CGAL::Polygon_mesh_processing::connected_components(g, boost::make_assoc_property_map(map), CGAL::Polygon_mesh_processing::parameters::all_default());
  Adapter fg(g, boost::make_assoc_property_map(map), 0);

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
void test_read(const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::face_descriptor g_face_descriptor;
  typedef boost::associative_property_map< std::map< g_face_descriptor, std::size_t > >FCMap;
  typedef CGAL::Connected_components_graph<Graph, FCMap> Adapter;
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  std::map<g_face_descriptor, std::size_t> map;
  CGAL::Polygon_mesh_processing::connected_components(g, boost::make_assoc_property_map(map), CGAL::Polygon_mesh_processing::parameters::all_default());
  Adapter fg(g, boost::make_assoc_property_map(map), 0);
  assert(CGAL::is_valid(fg));
}

template <typename Graph>
void
test(const std::vector<Graph>& graphs)
{
  BOOST_FOREACH(Graph p, graphs){
    test_read(p);
    test_vertex_iterators(p);
    test_out_edges(p);
    test_in_edges(p);
    test_in_out_edges(p);
    test_edge_find(p);
    test_faces(p);
    test_edge_iterators(p);
    test_halfedge_around_face_iterator(p);
    test_halfedge_around_vertex_iterator(p);
  }
}


typedef SM::Point Point_3;

template<typename VertexPointPMap>
struct Constraint : public boost::put_get_helper<bool,Constraint<VertexPointPMap> >
{
  typedef typename boost::graph_traits<SM>::edge_descriptor edge_descriptor;
  typedef boost::readable_property_map_tag      category;
  typedef bool                                  value_type;
  typedef bool                                  reference;
  typedef edge_descriptor                       key_type;


  Constraint()
    :g_(NULL)
  {}

  Constraint(SM& g, VertexPointPMap vpp)
    : g_(&g), vppmap(vpp)
  {}

  bool operator[](edge_descriptor e) const
  {
    const SM& g = *g_;
    if(
       (boost::get(vppmap, target(e, g)) == Point_3(1,1,1) ||
        boost::get(vppmap, source(e, g)) == Point_3(1,1,1)) &&
       (boost::get(vppmap, target(e, g)) == Point_3(0,0,0) ||
        boost::get(vppmap, source(e, g)) == Point_3(0,0,0))
       ||
       (boost::get(vppmap, target(e, g)) == Point_3(1,1,1) ||
        boost::get(vppmap, source(e, g)) == Point_3(1,1,1)) &&
       (boost::get(vppmap, target(e, g)) == Point_3(0,0,1) ||
        boost::get(vppmap, source(e, g)) == Point_3(0,0,1))
       ||
       (boost::get(vppmap, target(e, g)) == Point_3(0,0,1) ||
        boost::get(vppmap, source(e, g)) == Point_3(0,0,1)) &&
       (boost::get(vppmap, target(e, g)) == Point_3(0,0,0) ||
        boost::get(vppmap, source(e, g)) == Point_3(0,0,0))
      ||
      (boost::get(vppmap, target(e, g)) == Point_3(1,0,1) ||
       boost::get(vppmap, source(e, g)) == Point_3(1,0,1)) &&
      (boost::get(vppmap, target(e, g)) == Point_3(0,0,0) ||
       boost::get(vppmap, source(e, g)) == Point_3(0,0,0))
      ||
      (boost::get(vppmap, target(e, g)) == Point_3(1,1,1) ||
       boost::get(vppmap, source(e, g)) == Point_3(1,1,1)) &&
      (boost::get(vppmap, target(e, g)) == Point_3(1,0,1) ||
       boost::get(vppmap, source(e, g)) == Point_3(1,0,1))
       )
      return true;
    else
      return false;
  }

  const SM* g_;
  VertexPointPMap vppmap;
};


int
main()
{
  typedef CGAL::Connected_components_graph<SM, SM::Property_map<boost::graph_traits<SM>::face_descriptor , std::size_t> > Adapter;
 // test(sm_data());
  //Make a tetrahedron and test the adapter for a patch that only contains 2 faces
  SM* sm = new SM();
  CGAL::make_tetrahedron(
        Point_3(1,1,1),
        Point_3(0,0,0),
        Point_3(0,0,1),
        Point_3(1,0,1),
        *sm);
  SM::Property_map<boost::graph_traits<SM>::face_descriptor , std::size_t> fccmap =
      sm->add_property_map<boost::graph_traits<SM>::face_descriptor, std::size_t>("f:CC").first;
  CGAL::Properties::Property_map<boost::graph_traits<SM>::vertex_descriptor, SM::Point> positions =
      sm->points();
  CGAL::Polygon_mesh_processing::connected_components(*sm, fccmap, CGAL::Polygon_mesh_processing::parameters::
                                                      edge_is_constrained_map(Constraint<CGAL::Properties::Property_map<boost::graph_traits<SM>::vertex_descriptor,
                                                                              SM::Point> >(*sm, positions)));
  boost::unordered_set<long unsigned int> pids;
  pids.insert(0);
  pids.insert(2);
  Adapter fga(*sm, fccmap, pids);
  CGAL_GRAPH_TRAITS_MEMBERS(Adapter);
  //check that there is the right number of simplices in fga
  CGAL_assertion(CGAL::is_valid(fga));
  CGAL_assertion(std::distance(faces(fga).first,faces(fga).second) == 2);
  CGAL_assertion(std::distance(edges(fga).first,edges(fga).second) == 5);
  CGAL_assertion(std::distance(halfedges(fga).first,halfedges(fga).second) == 10);
  CGAL_assertion(std::distance(vertices(fga).first,vertices(fga).second) == 4);
  halfedge_descriptor h = halfedge(*faces(fga).first, fga);
  vertex_descriptor v = source(h, fga);
  //check that next() works inside the patch
  CGAL_assertion(
        next(next(next(h, fga), fga), fga) == h
        );
  //check that next() works on bordure of the patch
  h = opposite(h, fga);
  CGAL_assertion(
        next(next(next(next(h, fga), fga), fga), fga) == h
        );
  //check that prev() works inside the patch
   h = halfedge(*faces(fga).first, fga);
  CGAL_assertion(
        prev(prev(prev(h, fga), fga), fga) == h
        );
  //check that prev() works on bordure of the patch
  h = opposite(h, fga);
  CGAL_assertion(
        prev(prev(prev(prev(h, fga), fga), fga), fga) == h
        );
  //check degree
  CGAL_assertion(degree(v, fga) == 3);
  //check in_edges and out_edges
  CGAL_assertion(std::distance(in_edges(v, fga).first ,in_edges(v, fga).second) == 3 );
  CGAL_assertion(std::distance(out_edges(v, fga).first ,out_edges(v, fga).second) == 3 );

  return 0;
}
