#include "test_Prefix.h"

#include <boost/unordered_set.hpp>

template<typename G,
         typename ForwardRange,
         typename IndexPropertyMap>
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

template <typename Graph>
void index_uniqueness(const Graph& g)
{
  std::cout << "Graph has:" << std::endl
            << "\t" << num_vertices(g) << " vertices" << std::endl
            << "\t" << num_halfedges(g) << " halfedges" << std::endl
            << "\t" << num_edges(g) << " edges" << std::endl
            << "\t" << num_faces(g) << " faces" << std::endl;

  index_uniqueness(g, vertices(g),  get(boost::vertex_index, g));
  index_uniqueness(g, halfedges(g), get(boost::halfedge_index, g));
  index_uniqueness(g, edges(g) ,    get(boost::edge_index, g));
  index_uniqueness(g, faces(g),     get(boost::face_index, g));

  index_uniqueness(g, vertices(g),  CGAL::get_initialized_vertex_index_map(g));
  index_uniqueness(g, halfedges(g), CGAL::get_initialized_halfedge_index_map(g));
  index_uniqueness(g, edges(g) ,    CGAL::get_initialized_edge_index_map(g));
  index_uniqueness(g, faces(g),     CGAL::get_initialized_face_index_map(g));
}

void index_uniqueness_poly(const Polyhedron& g)
{
  index_uniqueness(g);

  index_uniqueness(g, edges(g) ,    get(boost::edge_external_index, g));
  index_uniqueness(g, vertices(g),  get(boost::vertex_external_index, g));
  index_uniqueness(g, faces(g),     get(boost::face_external_index, g));
  index_uniqueness(g, halfedges(g), get(boost::halfedge_external_index, g));
}

// @todo test non-const
int main(int, char**)
{
  std::cout << "testing Polyhedron\n";
  std::vector<Polyhedron> polys = poly_data(); // @tmp
  for(const Polyhedron& p : polys)
    index_uniqueness_poly(p);

  std::cout << "testing Linear_cell_complex\n";
  std::vector<LCC> lccs = lcc_data();
  for(const LCC& p : lccs)
    index_uniqueness(p);

  std::cout << "testing Surface_mesh\n";
  std::vector<SM> sms = sm_data();
  for(const SM& p : sms)
    index_uniqueness(p);

#if defined(CGAL_USE_OPENMESH)
  std::cout << "testing OpenMesh\n";
  std::vector<OMesh> omeshs = omesh_data();
  for(const OMesh& p : omeshs)
    index_uniqueness(p);
#endif

  std::cout << "testing Triangulations\n";
  index_uniqueness(t2_data());
  index_uniqueness(dt2_data());
  index_uniqueness(rt2_data());
  index_uniqueness(ct2_data());
  index_uniqueness(cdt2_data());
  index_uniqueness(cdtp2_data());
  index_uniqueness(t2h_data());

  std::cerr << "done\n";
  return 0;
}
