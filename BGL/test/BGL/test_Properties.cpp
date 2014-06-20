
#define BOOST_TEST_MAIN 1
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include "test_Prefix.h"
#include <boost/unordered_set.hpp>

template< typename G,
          typename ForwardRange,
          typename IndexPropertyMap
          >
void index_uniqueness(const G&, 
                      ForwardRange range,
                      IndexPropertyMap pm)
{
  typename boost::range_iterator<ForwardRange>::type 
    begin = boost::begin(range), 
    begin2 = boost::begin(range),
    end = boost::end(range);

  typedef boost::unordered_set<typename IndexPropertyMap::value_type> id_map;
  typedef std::pair<typename id_map::iterator, bool> resultp;
  id_map m;

  while(begin != end) {
    resultp r = m.insert(get(pm, *begin));
    ++begin;
    BOOST_CHECK(r.second);
  }

  BOOST_CHECK_EQUAL(std::distance(begin2, end), m.size());
}


void index_uniqueness_poly(const Polyhedron& g)
{
  index_uniqueness(g, edges(g) ,    get(boost::edge_index, g));
  index_uniqueness(g, vertices(g),  get(boost::vertex_index, g));
  index_uniqueness(g, faces(g),     get(boost::face_index, g));
  index_uniqueness(g, halfedges(g), get(boost::halfedge_index, g));

  index_uniqueness(g, edges(g) ,    get(boost::edge_external_index, g));
  index_uniqueness(g, vertices(g),  get(boost::vertex_external_index, g));
  index_uniqueness(g, faces(g),     get(boost::face_external_index, g));
  index_uniqueness(g, halfedges(g), get(boost::halfedge_external_index, g));
}
#if defined(CGAL_USE_SURFACE_MESH)
void index_uniqueness_sm(const SM& g)
{
  index_uniqueness(g, edges(g) ,    get(boost::edge_index, g));
  index_uniqueness(g, vertices(g),  get(boost::vertex_index, g));
  index_uniqueness(g, faces(g),     get(boost::face_index, g));
  index_uniqueness(g, halfedges(g), get(boost::halfedge_index, g));
}
#endif

#if defined(CGAL_USE_OPENMESH)
void index_uniqueness_omesh(const OMesh& g)
{
  index_uniqueness(g, edges(g) ,    get(boost::edge_index, g));
  index_uniqueness(g, vertices(g),  get(boost::vertex_index, g));
  index_uniqueness(g, faces(g),     get(boost::face_index, g));
  index_uniqueness(g, halfedges(g), get(boost::halfedge_index, g));
}
#endif

using namespace boost::unit_test;

test_suite*
init_unit_test_suite( int, char** const )
{
  std::vector<Polyhedron> polys = poly_data();


  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE(&index_uniqueness_poly, polys.begin(), polys.end() ) );

#if defined(CGAL_USE_SURFACE_MESH)
  std::vector<SM> sms = sm_data();

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE(&index_uniqueness_sm, sms.begin(), sms.end() ) );
#endif 

#if defined(CGAL_USE_OPENMESH)
  std::vector<OMesh> omeshs = omesh_data();
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE(&index_uniqueness_omesh, omeshs.begin(), omeshs.end() ) );
#endif

  return 0;
}


// trick cgal_test_with_cmake into adding this file to the test-suite
// int main()
// {
//   return 0;
// }
