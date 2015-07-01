

#include "test_Prefix.h"
#include <boost/numeric/conversion/cast.hpp>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>
#include <CGAL/use.h>

typedef boost::unordered_set<std::size_t> id_map;

template <typename Graph>
void test_isolated_vertex(const Graph& g)
{
  std::cerr << typeid(g).name() << std::endl;
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
      halfedge_around_target_iterator step = boost::next(havit);
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
  assert(boost::numeric_cast<edges_size_type>(std::distance(eb, ee)) == num_edges(g));

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
void test_vertex_iterators(G& g)
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
  assert(CGAL::is_valid(g));
}


template <typename Graph>
void
test(const std::vector<Graph>& graphs)
{
  BOOST_FOREACH(Graph p, graphs){
    test_edge_iterators(p);
    test_read(p);
    test_vertex_iterators(p);
    test_out_edges(p);
    test_in_edges(p);
    test_in_out_edges(p);
    test_edge_find(p);
    test_faces(p);
    test_halfedge_around_vertex_iterator(p);
    test_halfedge_around_face_iterator(p);
    test_isolated_vertex(p);
  }
}

int
main()
{
  test(poly_data());

#if defined(CGAL_USE_SURFACE_MESH)
  test(sm_data());
#endif

#if defined(CGAL_USE_OPENMESH)
  test(omesh_data());
#endif

  std::cerr << "done" << std::endl;
  return 0;
}
