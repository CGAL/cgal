
#include "test_Prefix.h"
#include <boost/unordered_set.hpp>
#include <boost/foreach.hpp>

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
    assert(r.second);
  }

  assert(std::distance(begin2, end) == static_cast<std::ptrdiff_t>(m.size()));
}


void index_uniqueness_poly(const Polyhedron& g)
{
  std::cerr << "testing Polyhedron\n";
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
  std::cerr << "testing Surface_mesh\n";
  index_uniqueness(g, edges(g) ,    get(boost::edge_index, g));
  index_uniqueness(g, vertices(g),  get(boost::vertex_index, g));
  index_uniqueness(g, faces(g),     get(boost::face_index, g));
  index_uniqueness(g, halfedges(g), get(boost::halfedge_index, g));
}
#endif

#if defined(CGAL_USE_OPENMESH)
void index_uniqueness_omesh(const OMesh& g)
{
  std::cerr << "testing OpenMesh\n";
  index_uniqueness(g, edges(g) ,    get(boost::edge_index, g));
  index_uniqueness(g, vertices(g),  get(boost::vertex_index, g));
  index_uniqueness(g, faces(g),     get(boost::face_index, g));
  index_uniqueness(g, halfedges(g), get(boost::halfedge_index, g));
}
#endif


int
main()
{
  std::vector<Polyhedron> polys = poly_data();

  BOOST_FOREACH(Polyhedron p, polys){
    index_uniqueness_poly(p);
  }



#if defined(CGAL_USE_SURFACE_MESH)
  std::vector<SM> sms = sm_data();

  BOOST_FOREACH(SM p, sms){
    index_uniqueness_sm(p);
  }
#endif 

#if defined(CGAL_USE_OPENMESH)
  std::vector<OMesh> omeshs = omesh_data();
  BOOST_FOREACH(OMesh p, omeshs){
    index_uniqueness_omesh(p);
  }
#endif

  std::cerr << "done\n";
  return 0;
}
