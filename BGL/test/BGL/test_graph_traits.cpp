
#define BOOST_TEST_MAIN 1

#include <boost/test/parameterized_test.hpp>
#include <boost/test/test_case_template.hpp>

#include "test_Prefix.h"

#include <boost/unordered_set.hpp>

typedef boost::unordered_set<std::size_t> id_map;

template <typename Graph>
void test_halfedge_around_vertex_iterator(const Graph& g)
{
  CGAL_GRAPH_TRAITS_MEMBERS(Graph);
  vertex_iterator vit, vend;
  for(boost::tie(vit, vend) = vertices(g); vit != vend; ++vit) {
    halfedge_around_target_iterator havit, havend;
    for(boost::tie(havit, havend) = CGAL::halfedges_around_target(halfedge(*vit, g), g); 
        havit != havend; ++havit) {
      BOOST_CHECK(target(*havit, g) == *vit);

      // check if we are really moving clockwise
      halfedge_around_target_iterator step = boost::next(havit);
      if(step != havend) {
        halfedge_descriptor stepd = *step;
        BOOST_CHECK(stepd == opposite(next(*havit, g), g));
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
    BOOST_CHECK_NE(std::distance(hafit, hafend), 0);
    for(boost::tie(hafit, hafend) = CGAL::halfedges_around_face(halfedge(*fit, g), g); hafit != hafend; ++hafit) {
      BOOST_CHECK(face(*hafit, g) == *fit);
    }
  }
}

template<typename G>
void test_edge_iterators(const G& g)
{
  typedef boost::graph_traits< G > Traits;
  typedef typename Traits::edge_descriptor edge_descriptor;
  typedef typename Traits::edge_iterator edge_iterator;

  BOOST_CHECK_EQUAL(g.size_of_halfedges() / 2, num_edges(g));
  
  // do we iterate as many as that?
  edge_iterator eb, ee;
  boost::tie(eb, ee) = edges(g);
  BOOST_CHECK(std::distance(eb, ee) == static_cast<std::ptrdiff_t>(g.size_of_halfedges() / 2));

  id_map ids;
  unsigned int count = 0;
  for(boost::tie(eb, ee) = edges(g); eb != ee; ++eb) {
    edge_descriptor e = *eb;
    std::pair<id_map::iterator, bool> r = ids.insert(get(boost::edge_index, g, e));
    // unique?
    BOOST_CHECK(r.second);
    ++count;
  }
  BOOST_CHECK_EQUAL(count, num_edges(g));
}

template<typename G>
void test_vertex_iterators(G& g)
{
  typedef boost::graph_traits< G > Traits;
  typedef typename Traits::vertex_iterator vertex_iterator;

  std::size_t count = 0;
  for (typename G::Vertex_iterator it = g.vertices_begin(); it != g.vertices_end(); ++it)
    ++count;

  BOOST_CHECK_EQUAL(count, num_vertices(g));

  // check that the iterators reach uniques
  id_map ids;
  vertex_iterator vb, ve;
  count = 0;
  for(boost::tie(vb, ve) = vertices(g); vb != ve; ++vb) {
    std::pair<id_map::iterator, bool> r = ids.insert(get(boost::vertex_index, g, *vb));
    BOOST_CHECK(r.second);
    ++count;
  }
  BOOST_CHECK_EQUAL(count, num_vertices(g));
}

template<typename G>
void test_out_edges(const G& g) 
{
  typedef boost::graph_traits< G > Traits;
  typedef typename Traits::vertex_iterator vertex_iterator;
  typedef typename Traits::in_edge_iterator in_edge_iterator;
  typedef typename Traits::out_edge_iterator out_edge_iterator;
  typedef typename Traits::vertex_descriptor vertex_descriptor;
  typedef typename Traits::edge_descriptor edge_descriptor;

  vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = vertices(g); vb != ve; ++vb) {
    id_map v_ids;

    vertex_descriptor around = *vb;
    out_edge_iterator oeb, oee;
    for(boost::tie(oeb, oee) = out_edges(*vb, g); oeb != oee; ++oeb) {
      vertex_descriptor t = target(*oeb, g);
      vertex_descriptor s = source(*oeb, g);
      BOOST_CHECK(s != t);
      BOOST_CHECK(s == around);
      BOOST_CHECK(t != around);
      std::pair<id_map::iterator, bool> r =
        v_ids.insert(get(boost::vertex_index, g, target(*oeb, g)));
      BOOST_CHECK(r.second);
    }
  }
}

template<typename G>
void test_in_edges(const G& g) 
{
  typedef boost::graph_traits< G > Traits;
  typedef typename Traits::vertex_iterator vertex_iterator;
  typedef typename Traits::in_edge_iterator in_edge_iterator;
  typedef typename Traits::out_edge_iterator out_edge_iterator;
  typedef typename Traits::vertex_descriptor vertex_descriptor;
  typedef typename Traits::edge_descriptor edge_descriptor;

  vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = vertices(g); vb != ve; ++vb) {
    id_map v_ids;
    vertex_descriptor around = *vb;
    in_edge_iterator ieb, iee;
    for(boost::tie(ieb, iee) = in_edges(*vb, g); ieb != iee; ++ieb) {
      vertex_descriptor t = target(*ieb, g);
      vertex_descriptor s = source(*ieb, g);
      BOOST_CHECK(t == around);
      BOOST_CHECK(s != around);
      std::pair<id_map::iterator, bool> r =
        v_ids.insert(get(boost::vertex_index, g, source(*ieb, g)));
      BOOST_CHECK(r.second);
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
  typedef typename Traits::edge_descriptor edge_descriptor;

  // check that the sets of in out edges are the same
  vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = vertices(g); vb != ve; ++vb) {
    id_map v_ids;
    std::vector<vertex_descriptor> in, out;
    in_edge_iterator ieb, iee;
    for(boost::tie(ieb, iee) = in_edges(*vb, g); ieb != iee; ++ieb) {
      std::pair<id_map::iterator, bool> r =
        v_ids.insert(get(boost::vertex_index, g, source(*ieb, g)));
      BOOST_CHECK(r.second);
      in.push_back(source(*ieb, g));
    }
    out_edge_iterator oeb, oee;
    for(boost::tie(oeb, oee) = out_edges(*vb, g); oeb != oee; ++oeb) {
      std::pair<id_map::iterator, bool> r =
        v_ids.insert(get(boost::vertex_index, g, target(*oeb, g)));
      // insertion must fail
      BOOST_CHECK(!r.second);
      out.push_back(target(*oeb, g));
    }
    // did we walk the vertices in the same order?
    BOOST_CHECK(in.size() == out.size());
    BOOST_CHECK(std::equal(in.begin(), in.end(), out.begin()));
    BOOST_CHECK(in.size() == in_degree(*vb, g));
    BOOST_CHECK(out.size() == out_degree(*vb, g));
    BOOST_CHECK(in.size() == degree(*vb, g));
    BOOST_CHECK(degree(*vb, g) == in_degree(*vb, g));
    BOOST_CHECK(degree(*vb, g) == out_degree(*vb, g));
  }
}

// check that every edge can be found through edge(u, v, g)
template<typename G>
void test_edge_find(const G& g)
{
  typedef boost::graph_traits<G>    Traits;
  typedef typename Traits::vertex_iterator   vertex_iterator;
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
    BOOST_CHECK(found.second);
    BOOST_CHECK(found2.second);
    BOOST_CHECK(found.first == *eb);
    BOOST_CHECK(found2.first == *eb);
  }
}

template<typename G>
void test_faces(const G& g)
{
  typedef boost::graph_traits<G>                         Traits;
  typedef typename Traits::vertex_iterator               vertex_iterator;
  typedef typename Traits::edge_iterator                 edge_iterator;
  typedef typename Traits::face_iterator                 face_iterator;
  typedef typename Traits::vertex_descriptor             vertex_descriptor;
  typedef typename Traits::halfedge_descriptor           halfedge_descriptor;
  typedef CGAL::Halfedge_around_face_iterator<G>         halfedge_around_face_iterator;

  unsigned int count = 0;
  face_iterator fb, fe;
  for(boost::tie(fb, fe) = faces(g); fb != fe; ++fb) {
    ++count;
    // reverse look-up
    halfedge_descriptor assoc = halfedge(*fb, g);
    BOOST_CHECK(face(assoc, g) == *fb);
    // check the enclosure
    halfedge_around_face_iterator encb, ence;
    for(boost::tie(encb, ence) = CGAL::halfedges_around_face(halfedge(*fb, g), g); encb != ence; ++encb) {
      BOOST_CHECK(face(*encb, g) == *fb);
    }
  }
  BOOST_CHECK_EQUAL(count, num_faces(g));
}

template<typename G>
void test_read(const G& g)
{
  BOOST_CHECK(CGAL::is_valid(g));
}

using namespace boost::unit_test;

test_suite*
init_unit_test_suite( int , char** const)
{
  std::vector<Polyhedron> polys = poly_data();

  framework::master_test_suite().p_name.value = "test_graph_traits test suite";

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_edge_iterators<Polyhedron>, polys.begin(), polys.end() ) );

  // Polyhedron
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_read<Polyhedron>, polys.begin(), polys.end() ) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_vertex_iterators<Polyhedron>, polys.begin(), polys.end() ) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_out_edges<Polyhedron>, polys.begin(), polys.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_in_edges<Polyhedron>, polys.begin(), polys.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_in_out_edges<Polyhedron>, polys.begin(), polys.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_edge_find<Polyhedron>, polys.begin(), polys.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_faces<Polyhedron>, polys.begin(), polys.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_halfedge_around_vertex_iterator<Polyhedron>, polys.begin(), polys.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_halfedge_around_face_iterator<Polyhedron>, polys.begin(), polys.end()) );
#if defined(CGAL_USE_SURFACE_MESH)
  // Surface_mesh
  std::vector<SM> sms = sm_data();

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_vertex_iterators<SM>, sms.begin(), sms.end() ) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_out_edges<SM>, sms.begin(), sms.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_in_edges<SM>, sms.begin(), sms.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_in_out_edges<SM>, sms.begin(), sms.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_edge_find<SM>, sms.begin(), sms.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_faces<SM>, sms.begin(), sms.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_read<SM>, sms.begin(), sms.end() ) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_halfedge_around_vertex_iterator<SM>, sms.begin(), sms.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_halfedge_around_face_iterator<SM>, sms.begin(), sms.end()) );
#endif 

#if defined(CGAL_USE_OPENMESH)
  std::vector<OMesh> omeshs = omesh_data();
  // framework::master_test_suite().
  //   add( BOOST_PARAM_TEST_CASE( &test_vertex_iterators<OMesh>, omeshs.begin(), omeshs.end() ) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_out_edges<OMesh>, omeshs.begin(), omeshs.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_in_edges<OMesh>, omeshs.begin(), omeshs.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_in_out_edges<OMesh>, omeshs.begin(), omeshs.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_edge_find<OMesh>, omeshs.begin(), omeshs.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_faces<OMesh>, omeshs.begin(), omeshs.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_read<OMesh>, omeshs.begin(), omeshs.end() ) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_halfedge_around_vertex_iterator<OMesh>, omeshs.begin(), omeshs.end()) );
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &test_halfedge_around_face_iterator<OMesh>, omeshs.begin(), omeshs.end()) );
#endif

  return 0;
}


// trick cgal_test_with_cmake into adding this file to the test-suite
// int main()
// {
//   return 0;
// }

