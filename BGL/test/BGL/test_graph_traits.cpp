#include "test_Prefix.h"

#include <CGAL/use.h>

#include <unordered_set>

typedef std::unordered_set<std::size_t>                        id_map;

template <typename Graph>
void test_isolated_vertex()
{
  Graph G;
  typedef boost::graph_traits< Graph > Traits;
  typedef typename Traits::vertex_descriptor vertex_descriptor;
  typedef typename Traits::halfedge_descriptor halfedge_descriptor;
  vertex_descriptor v = add_vertex(G);
  // the connectivity of v may be anything
  set_halfedge(v, Traits::null_halfedge(), G);
  halfedge_descriptor h = halfedge(v,G);
  CGAL_USE(h);
}


template <typename Graph>
void test_halfedge_around_vertex_iterator(const Graph& g)
{
  CGAL_GRAPH_TRAITS_MEMBERS(Graph);
  vertex_iterator vit, vend;
  for(boost::tie(vit, vend) = vertices(g); vit != vend; ++vit) {
    halfedge_around_target_iterator havit, havend;
    for(boost::tie(havit, havend) = CGAL::halfedges_around_target(halfedge(*vit, g), g);
        havit != havend; ++havit) {
      assert(target(*havit, g) == *vit);

      // check if we are really moving clockwise
      halfedge_around_target_iterator step = std::next(havit);
      if(step != havend) {
        halfedge_descriptor stepd = *step;
        assert(stepd == opposite(next(*havit, g), g));
      }
    }
  }
}

template <typename Graph>
void test_halfedge_around_face_iterator(const Graph& g)
{
  CGAL_GRAPH_TRAITS_MEMBERS(Graph);
  face_iterator fit, fend;
  for(boost::tie(fit, fend) = faces(g); fit != fend; ++fit) {
    halfedge_around_face_iterator hafit, hafend;
    boost::tie(hafit, hafend) = CGAL::halfedges_around_face(halfedge(*fit, g), g);
    assert(std::distance(hafit, hafend) != 0);
    for(boost::tie(hafit, hafend) = CGAL::halfedges_around_face(halfedge(*fit, g), g); hafit != hafend; ++hafit) {
      assert(face(*hafit, g) == *fit);
    }
  }
}

template<typename G>
void test_halfedge_iterators(const G& g)
{
  typedef boost::graph_traits< G > Traits;
  typedef typename Traits::halfedge_iterator halfedge_iterator;
  typedef typename Traits::halfedges_size_type halfedges_size_type;

  // do we iterate as many as that?
  halfedge_iterator hb, he;
  boost::tie(hb, he) = halfedges(g);
  assert(static_cast<halfedges_size_type>(std::distance(hb, he)) == num_halfedges(g));

  id_map ids;
  unsigned int count = 0;
  for(boost::tie(hb, he) = halfedges(g); hb != he; ++hb) {
    std::pair<id_map::iterator, bool> r = ids.insert(get(boost::halfedge_index, g, *hb));
    // unique?
    assert(r.second);
    ++count;
  }
  assert(count == num_halfedges(g));
}

template<typename G>
void test_edge_iterators(const G& g)
{
  typedef boost::graph_traits< G > Traits;
  typedef typename Traits::edge_descriptor edge_descriptor;
  typedef typename Traits::edge_iterator edge_iterator;
  typedef typename Traits::edges_size_type edges_size_type;

  // assert(g.size_of_halfedges() / 2 == num_edges(g));

  // do we iterate as many as that?
  edge_iterator eb, ee;
  boost::tie(eb, ee) = edges(g);
  assert(static_cast<edges_size_type>(std::distance(eb, ee)) == num_edges(g));

  id_map ids;
  unsigned int count = 0;
  for(boost::tie(eb, ee) = edges(g); eb != ee; ++eb) {
    edge_descriptor e = *eb;
    std::pair<id_map::iterator, bool> r = ids.insert(get(boost::edge_index, g, e));
    // unique?
    assert(r.second);
    ++count;
  }
  assert(count == num_edges(g));
}

template<typename G>
void test_vertex_iterators(const G& g)
{
  typedef boost::graph_traits< G > Traits;
  typedef typename Traits::vertex_iterator vertex_iterator;

  vertex_iterator vb, ve;
  std::size_t count = 0;
  for(boost::tie(vb, ve) = vertices(g); vb != ve; ++vb){
    ++count;
  }

  assert(count == num_vertices(g));

  // check that the iterators reach uniques
  id_map ids;

  count = 0;
  for(boost::tie(vb, ve) = vertices(g); vb != ve; ++vb) {
    std::pair<id_map::iterator, bool> r = ids.insert(get(boost::vertex_index, g, *vb));
    assert(r.second);
    ++count;
  }
  assert(count == num_vertices(g));
}

template<typename G>
void test_out_edges(const G& g)
{
  typedef boost::graph_traits< G > Traits;
  typedef typename Traits::vertex_iterator vertex_iterator;
  typedef typename Traits::out_edge_iterator out_edge_iterator;
  typedef typename Traits::vertex_descriptor vertex_descriptor;

  vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = vertices(g); vb != ve; ++vb) {
    id_map v_ids;

    vertex_descriptor around = *vb;
    out_edge_iterator oeb, oee;
    for(boost::tie(oeb, oee) = out_edges(*vb, g); oeb != oee; ++oeb) {
      vertex_descriptor t = target(*oeb, g);
      vertex_descriptor s = source(*oeb, g);
      assert(s != t);
      assert(s == around);
      assert(t != around);
      std::pair<id_map::iterator, bool> r =
        v_ids.insert(get(boost::vertex_index, g, target(*oeb, g)));
      assert(r.second);
    }
  }
}

template<typename G>
void test_in_edges(const G& g)
{
  typedef boost::graph_traits< G > Traits;
  typedef typename Traits::vertex_iterator vertex_iterator;
  typedef typename Traits::in_edge_iterator in_edge_iterator;
  typedef typename Traits::vertex_descriptor vertex_descriptor;

  vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = vertices(g); vb != ve; ++vb) {
    id_map v_ids;
    vertex_descriptor around = *vb;
    in_edge_iterator ieb, iee;
    for(boost::tie(ieb, iee) = in_edges(*vb, g); ieb != iee; ++ieb) {
      vertex_descriptor t = target(*ieb, g);
      vertex_descriptor s = source(*ieb, g);
      assert(t == around);
      assert(s != around);
      std::pair<id_map::iterator, bool> r =
        v_ids.insert(get(boost::vertex_index, g, source(*ieb, g)));
      assert(r.second);
    }
  }
}

template<typename G>
void test_in_out_edges(const G& g)
{
  typedef boost::graph_traits< G > Traits;
  typedef typename Traits::vertex_iterator vertex_iterator;
  typedef typename Traits::in_edge_iterator in_edge_iterator;
  typedef typename Traits::out_edge_iterator out_edge_iterator;
  typedef typename Traits::vertex_descriptor vertex_descriptor;

  // check that the sets of in out edges are the same
  vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = vertices(g); vb != ve; ++vb) {
    id_map v_ids;
    std::vector<vertex_descriptor> in, out;
    in_edge_iterator ieb, iee;
    for(boost::tie(ieb, iee) = in_edges(*vb, g); ieb != iee; ++ieb) {
      std::pair<id_map::iterator, bool> r =
        v_ids.insert(get(boost::vertex_index, g, source(*ieb, g)));
      assert(r.second);
      in.push_back(source(*ieb, g));
    }
    out_edge_iterator oeb, oee;
    for(boost::tie(oeb, oee) = out_edges(*vb, g); oeb != oee; ++oeb) {
      std::pair<id_map::iterator, bool> r =
        v_ids.insert(get(boost::vertex_index, g, target(*oeb, g)));
      // insertion must fail
      assert(!r.second);
      out.push_back(target(*oeb, g));
    }
    // did we walk the vertices in the same order?
    assert(in.size() == out.size());
    assert(std::equal(in.begin(), in.end(), out.begin()));
    assert(in.size() == in_degree(*vb, g));
    assert(out.size() == out_degree(*vb, g));
    assert(in.size() == degree(*vb, g));
    assert(degree(*vb, g) == in_degree(*vb, g));
    assert(degree(*vb, g) == out_degree(*vb, g));
  }
}

template<typename G>
void test_adjacent_vertices(const G& g)
{
  typedef boost::graph_traits< G > Traits;
  typedef typename Traits::vertex_descriptor vertex_descriptor;
  typedef typename Traits::edge_descriptor edge_descriptor;
  typedef typename Traits::in_edge_iterator in_edge_iterator;
  typedef typename Traits::out_edge_iterator out_edge_iterator;
  typedef typename Traits::adjacency_iterator adjacency_iterator;
  typedef std::pair<edge_descriptor, bool>   ret;

  vertex_descriptor v = *(vertices(g).begin());

  adjacency_iterator vb, ve;
  boost::tie(vb, ve) = adjacent_vertices(v, g);

  in_edge_iterator ieb, iee;
  boost::tie(ieb, iee) = in_edges(v, g);

  out_edge_iterator oeb, oee;
  boost::tie(oeb, oee) = out_edges(v, g);

  assert(std::distance(vb, ve) == std::distance(ieb, iee));
  assert(std::distance(vb, ve) == std::distance(oeb, oee));

  for(; vb != ve; ++vb)
  {
    vertex_descriptor s = *vb;
    assert(s != v);
    assert(s != Traits::null_vertex());
    ret found = edge(s, v, g);
    assert(found.second);
  }
}

// check that every edge can be found through edge(u, v, g)
template<typename G>
void test_edge_find(const G& g)
{
  typedef boost::graph_traits<G>    Traits;
  typedef typename Traits::edge_iterator     edge_iterator;
  typedef typename Traits::vertex_descriptor vertex_descriptor;
  typedef typename Traits::edge_descriptor   edge_descriptor;
  typedef std::pair<edge_descriptor, bool>   ret;

  edge_iterator eb, ee;
  for(boost::tie(eb, ee) = edges(g); eb != ee; ++eb) {
    vertex_descriptor s = source(*eb, g);
    vertex_descriptor t = target(*eb, g);
    ret found = edge(s, t, g);
    ret found2 = edge(t, s, g);
    assert(found.second);
    assert(found2.second);
    assert(found.first == *eb);
    assert(found2.first == *eb);
  }
}

template<typename G>
void test_faces(const G& g)
{
  typedef boost::graph_traits<G>                         Traits;
  typedef typename Traits::face_iterator                 face_iterator;
  typedef typename Traits::halfedge_descriptor           halfedge_descriptor;
  typedef CGAL::Halfedge_around_face_iterator<G>         halfedge_around_face_iterator;

  unsigned int count = 0;
  face_iterator fb, fe;
  for(boost::tie(fb, fe) = faces(g); fb != fe; ++fb) {
    ++count;
    // reverse look-up
    halfedge_descriptor assoc = halfedge(*fb, g);
    assert(face(assoc, g) == *fb);
    // check the enclosure
    halfedge_around_face_iterator encb, ence;
    for(boost::tie(encb, ence) = CGAL::halfedges_around_face(halfedge(*fb, g), g); encb != ence; ++encb) {
      assert(face(*encb, g) == *fb);
    }
  }
  assert(count == num_faces(g));
}

template<typename G>
void test_read(const G& g)
{
  assert(CGAL::is_valid_polygon_mesh(g));
}

template <typename Graph>
void test_const_graph(const Graph& g)
{
  test_vertex_iterators(g);
  test_halfedge_iterators(g);
  test_edge_iterators(g);
  test_read(g);
  test_out_edges(g);
  test_in_edges(g);
  test_in_out_edges(g);
//  test_adjacent_vertices(g);
  test_edge_find(g);
  test_faces(g);
  test_halfedge_around_vertex_iterator(g);
  test_halfedge_around_face_iterator(g);
}

template <typename Graph>
void test_graph_range(const std::vector<Graph>& graphs)
{
  for(const Graph& g : graphs)
  {
    test_const_graph(g);
    test_isolated_vertex<Graph>();
  }
}

int main()
{
  std::cout << "Test polyhedron data..." << std::endl;
  test_graph_range(poly_data());

  std::cout << "Test LCC data..." << std::endl;
  test_graph_range(lcc_data());

  std::cout << "Test Surface_mesh data..." << std::endl;
  test_graph_range(sm_data());

  std::cout << "Test T2 data..." << std::endl;
  test_const_graph(t2_data());
  test_const_graph(dt2_data());
  test_const_graph(rt2_data());
  test_const_graph(ct2_data());
  test_const_graph(cdt2_data());
  test_const_graph(cdtp2_data());
  test_const_graph(t2h_data());

#if defined(CGAL_USE_OPENMESH)
  std::cout << "Test OpenMesh data..." << std::endl;
  test_graph_range(omesh_data());
#endif

  std::cerr << "done" << std::endl;
  return 0;
}
