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
    assert(r.second);
  }

  assert(std::distance(begin2, end) == static_cast<std::ptrdiff_t>(m.size()));
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

void index_uniqueness_lcc(const LCC& g)
{
  index_uniqueness(g, edges(g) ,    get(boost::edge_index, g));
  index_uniqueness(g, vertices(g),  get(boost::vertex_index, g));
  index_uniqueness(g, faces(g),     get(boost::face_index, g));
  index_uniqueness(g, halfedges(g), get(boost::halfedge_index, g));
}

void index_uniqueness_sm(const SM& g)
{
  index_uniqueness(g, edges(g) ,    get(boost::edge_index, g));
  index_uniqueness(g, vertices(g),  get(boost::vertex_index, g));
  index_uniqueness(g, faces(g),     get(boost::face_index, g));
  index_uniqueness(g, halfedges(g), get(boost::halfedge_index, g));
}

#if defined(CGAL_USE_OPENMESH)
void index_uniqueness_omesh(const OMesh& g)
{
  index_uniqueness(g, edges(g) ,    get(boost::edge_index, g));
  index_uniqueness(g, vertices(g),  get(boost::vertex_index, g));
  index_uniqueness(g, faces(g),     get(boost::face_index, g));
  index_uniqueness(g, halfedges(g), get(boost::halfedge_index, g));
}
#endif

template <typename Triangulation>
void index_uniqueness_tr(const Triangulation& g)
{
  index_uniqueness(g, edges(g) ,    get(boost::edge_index, g));
  index_uniqueness(g, vertices(g),  get(boost::vertex_index, g));
  index_uniqueness(g, faces(g),     get(boost::face_index, g));
  index_uniqueness(g, halfedges(g), get(boost::halfedge_index, g));
}

int main()
{
  std::cout << "testing Polyhedron\n";
  std::vector<Polyhedron> polys = poly_data();
  for(Polyhedron p : polys)
    index_uniqueness_poly(p);

  std::cout << "testing Linear_cell_complex\n";
  std::vector<LCC> lccs = lcc_data();
  for(LCC p : lccs)
    index_uniqueness_lcc(p);

  std::cout << "testing Surface_mesh\n";
  std::vector<SM> sms = sm_data();
  for(SM p : sms)
    index_uniqueness_sm(p);

#if defined(CGAL_USE_OPENMESH)
  std::cout << "testing OpenMesh\n";
  std::vector<OMesh> omeshs = omesh_data();
  for(OMesh p : omeshs)
    index_uniqueness_omesh(p);
#endif

  std::cout << "testing Triangulations\n";
  index_uniqueness_tr(t2_data());
  index_uniqueness_tr(dt2_data());
  index_uniqueness_tr(rt2_data());
  index_uniqueness_tr(ct2_data());
  index_uniqueness_tr(cdt2_data());
  index_uniqueness_tr(cdtp2_data());
  index_uniqueness_tr(t2h_data());

  std::cerr << "done\n";
  return 0;
}
