#define BOOST_TEST_MODULE graph_traits test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_CASE( edges_test )
{
  edge_iterator eb, ee;
  vertex_iterator vb, ve;
  
  Surface_fixture f;
  boost::tie(eb, ee) = edges(f.m);
  boost::tie(vb, ve) = vertices(f.m);
  BOOST_CHECK(std::distance(eb, ee) == 7);
  BOOST_CHECK(std::distance(vb, ve) == 5);

  Cube_fixture cf;
  boost::tie(eb, ee) = edges(cf.m);
  boost::tie(vb, ve) = vertices(cf.m);
  BOOST_CHECK(std::distance(eb, ee) == 18);
  BOOST_CHECK(std::distance(vb, ve) == 8);
}

BOOST_AUTO_TEST_CASE( out_edges_test )
{
  Surface_fixture f;

  BOOST_CHECK(out_degree(f.u, f.m) == 2);
  BOOST_CHECK(out_degree(f.v, f.m) == 4);
  std::pair<out_edge_iterator, out_edge_iterator> out = out_edges(f.v, f.m);
  BOOST_CHECK(std::distance(out.first, out.second) == 4);
  out_edge_iterator it = out.first;
  while(it != out.second) {
    // the source should always be u
    BOOST_CHECK(source(*it, f.m) == f.v);
    // the target never
    BOOST_CHECK(target(*it, f.m) != f.v);
   ++it;
  }
}

BOOST_AUTO_TEST_CASE( in_edges_test )
{
  Surface_fixture f;
  BOOST_CHECK(in_degree(f.u, f.m) == 2);
  BOOST_CHECK(in_degree(f.x, f.m) == 3);
  BOOST_CHECK(in_degree(f.v, f.m) == 4);
  std::pair<in_edge_iterator, in_edge_iterator> in = in_edges(f.v, f.m);
  BOOST_CHECK(std::distance(in.first, in.second) == 4);

  in_edge_iterator it = in.first;
  while(it != in.second) {
    // the source should never be u
    BOOST_CHECK(source(*it, f.m) != f.v);
    // the target must always be u
    BOOST_CHECK(target(*it, f.m) == f.v);
    ++it;
  }
}

BOOST_AUTO_TEST_CASE( in_out_equality )
{
  // in and out degrees must be equal for each vertex
  Cube_fixture f;
  for(Sm::Vertex_iterator it = f.m.vertices_begin();
      it != f.m.vertices_end(); ++it) {
    BOOST_CHECK(in_degree(*it, f.m) == out_degree(*it, f.m));
  }
}

BOOST_AUTO_TEST_CASE( face_test )
{
  Surface_fixture f;
  std::pair<enclosure_iterator, enclosure_iterator> 
    enc = enclosure(f.f1, f.m);
  BOOST_CHECK(enc.first != enc.second);
  BOOST_CHECK(std::distance(enc.first, enc.second) == 3);
  enclosure_iterator begin = enc.first;
  while(begin != enc.second) 
  {
    BOOST_CHECK(face(*begin, f.m) == f.f1);
    ++begin;
  }
}

BOOST_AUTO_TEST_CASE( weight_map_test )
{
  Surface_fixture f;
  Cube_fixture c;
  
  CGAL::SM_edge_weight_pmap<K> wm1 = boost::get(boost::edge_weight, c.m);
  edge_iterator eb, ee;
  boost::test_tools::check_is_close_t check_close;
  for(boost::tie(eb, ee) = edges(c.m); eb != ee; ++eb) {
    BOOST_CHECK(
      check_close(wm1[*eb], 2.0, boost::test_tools::percent_tolerance_t<double>(0.00001))
      || check_close(wm1[*eb], 2.82843, boost::test_tools::percent_tolerance_t<double>(0.001)));
  }
}


BOOST_AUTO_TEST_CASE( vertices_test )
{
  Surface_fixture f;

  // boost::property_map<Sm,boost::vertex_index_t>::type vi_map = get(f.m, boost::vertex_index);

  // vertex_iterator b,e;

  // for(boost::tie(b,e) = vertices(f.m);
  //     b!= e;
  //     ++b){
  //   std::cout << boost::get(vi_map, *(b)) << std::endl;
  // }

  // boost::property_map<Sm,boost::edge_weight_t>::type ew_map = get(f.m, boost::edge_weight);

  // edge_iterator be, ee;

  // for(boost::tie(be,ee) = edges(f.m);
  //     be!= ee;
  //     ++be){
  //   std::cout << boost::get(ew_map, *(be)) << std::endl;
  // }

}
