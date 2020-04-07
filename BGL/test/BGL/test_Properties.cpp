#include "test_Prefix.h"

#include <CGAL/boost/graph/Euler_operations.h>

#include <boost/unordered_set.hpp>

// #define CGAL_TEST_PROPERTIES_DEBUG

namespace CGAL {

template <typename Key, typename Value, typename Container>
struct Non_mutable_property_map
{
  typedef Key                                     key_type;
  typedef Value                                   value_type;
  typedef value_type                              reference;
  typedef boost::readable_property_map_tag        category;

  Non_mutable_property_map(const Container& c) : m_c(c) { }

  friend reference get(const Non_mutable_property_map<Key, Value, Container>& pmap, key_type k)
  {
    return pmap.m_c.at(k);
  }

private:
  const Container& m_c;
};

template <typename Key, typename Value, typename Container>
struct RW_property_map
{
  typedef Key                                     key_type;
  typedef Value                                   value_type;
  typedef value_type&                             reference;
  typedef boost::read_write_property_map_tag      category;

  RW_property_map(Container& c) : m_c(c) { }

  friend void put(RW_property_map<Key, Value, Container>& pmap, const key_type& k, const value_type& val)
  {
    pmap.m_c[k] = val;
  }

  friend reference get(RW_property_map<Key, Value, Container>& pmap, const key_type& k)
  {
    return pmap.m_c[k];
  }

private:
  Container& m_c;
};

} // namespace CGAL

template<typename Graph,
         typename ForwardRange,
         typename IndexPropertyMap>
void test_uniqueness(const Graph&,
                     const ForwardRange& range,
                     IndexPropertyMap index_map)
{
#ifdef CGAL_TEST_PROPERTIES_DEBUG
  std::cout << std::endl
            << "Checking the uniqueness of the property map of type: "
            << typeid(IndexPropertyMap).name() << std::endl;
  std::cout << "Element type: " << typeid(typename boost::range_value<ForwardRange>::type).name() << std::endl;
#endif

  typename boost::range_iterator<ForwardRange>::type
    begin = boost::begin(range),
    begin2 = boost::begin(range),
    end = boost::end(range);

  typedef boost::unordered_set<typename IndexPropertyMap::value_type> id_map;
  typedef std::pair<typename id_map::iterator, bool> resultp;

  id_map m;
  while(begin != end)
  {
    resultp r = m.insert(get(index_map, *begin));
#ifdef CGAL_TEST_PROPERTIES_DEBUG
    std::cout << "id: " << get(index_map, *begin) << std::endl;
#endif
    ++begin;
    assert(r.second); // already seen that id
  }

  assert(std::distance(begin2, end) == static_cast<std::ptrdiff_t>(m.size()));
}

template<typename Graph,
         typename NamedParameters>
void test_vertex_index_map_uniqueness(const Graph& g,
                                      const NamedParameters& np)
{
  typedef typename CGAL::GetInitializedVertexIndexMap<Graph, NamedParameters>::type       VIM;
  typedef typename CGAL::GetInitializedVertexIndexMap<Graph, NamedParameters>::const_type CVIM;

  // in the case where the map is passed by NP, its type doesn't depend on whether the mesh is const or not
  static_assert((std::is_same<VIM, CVIM>::value), "VIM, CVIM must be the same type");

  VIM ivim = CGAL::get_initialized_vertex_index_map(g, np);

  return test_uniqueness(g, vertices(g), ivim);
}

template<typename Graph,
         typename NamedParameters>
void test_halfedge_index_map_uniqueness(const Graph& g,
                                        const NamedParameters& np)
{
  typedef typename CGAL::GetInitializedHalfedgeIndexMap<Graph, NamedParameters>::type       HIM;
  typedef typename CGAL::GetInitializedHalfedgeIndexMap<Graph, NamedParameters>::const_type CHIM;

  // in the case where the map is passed by NP, its type doesn't depend on whether the mesh is const or not
  static_assert((std::is_same<HIM, CHIM>::value), "HIM, CHIM must be the same type");

  HIM ihim = CGAL::get_initialized_halfedge_index_map(g, np);

  return test_uniqueness(g, halfedges(g), ihim);
}

template<typename Graph,
         typename NamedParameters>
void test_edge_index_map_uniqueness(const Graph& g,
                                    const NamedParameters& np)
{
  typedef typename CGAL::GetInitializedEdgeIndexMap<Graph, NamedParameters>::type       EIM;
  typedef typename CGAL::GetInitializedEdgeIndexMap<Graph, NamedParameters>::const_type CEIM;

  // in the case where the map is passed by NP, its type doesn't depend on whether the mesh is const or not
  static_assert((std::is_same<EIM, CEIM>::value), "EIM, CEIM must be the same type");

  EIM ieim = CGAL::get_initialized_edge_index_map(g, np);

  return test_uniqueness(g, edges(g), ieim);
}

template<typename Graph,
         typename NamedParameters>
void test_face_index_map_uniqueness(const Graph& g,
                                    const NamedParameters& np)
{
  typedef typename CGAL::GetInitializedFaceIndexMap<Graph, NamedParameters>::type       FIM;
  typedef typename CGAL::GetInitializedFaceIndexMap<Graph, NamedParameters>::const_type CFIM;

  // in the case where the map is passed by NP, its type doesn't depend on whether the mesh is const or not
  static_assert((std::is_same<FIM, CFIM>::value), "FIM, CFIM must be the same type");

  FIM ifim = CGAL::get_initialized_face_index_map(g, np);

  return test_uniqueness(g, faces(g), ifim);
}

////////////////////////////////////////// const ///////////////////////////////////////////////////

template <typename Graph>
void test_internal_index_maps_const(const Graph& g)
{
  test_uniqueness(g, vertices(g),  get(boost::vertex_index, g));
  test_uniqueness(g, halfedges(g), get(boost::halfedge_index, g));
  test_uniqueness(g, edges(g) ,    get(boost::edge_index, g));
  test_uniqueness(g, faces(g),     get(boost::face_index, g));
}

template <typename Graph>
void test_initialized_index_maps_const(const Graph& g)
{
  typedef typename CGAL::GetInitializedVertexIndexMap<Graph>::const_type      VIM;
  VIM ivim = CGAL::get_initialized_vertex_index_map(g);
  test_uniqueness(g, vertices(g), ivim);

  typedef typename CGAL::GetInitializedHalfedgeIndexMap<Graph>::const_type    HIM;
  HIM ihim = CGAL::get_initialized_halfedge_index_map(g);
  test_uniqueness(g, halfedges(g), ihim);

  typedef typename CGAL::GetInitializedEdgeIndexMap<Graph>::const_type        EIM;
  EIM ieim = CGAL::get_initialized_edge_index_map(g);
  test_uniqueness(g, edges(g), ieim);

  typedef typename CGAL::GetInitializedFaceIndexMap<Graph>::const_type        FIM;
  FIM ifim = CGAL::get_initialized_face_index_map(g);
  test_uniqueness(g, faces(g), ifim);

  // Passing an index map via NP
  typedef typename boost::graph_traits<Graph>::vertex_descriptor        vertex_descriptor;
  typedef std::map<vertex_descriptor, int>                              VertexIndexMap;
  typedef boost::associative_property_map<VertexIndexMap>               VertexIdPropertyMap; // lvalue_pmap

  int vi = static_cast<int>(num_vertices(g));
  VertexIndexMap vim;
  VertexIdPropertyMap external_vertex_index_map(vim);
  for(vertex_descriptor v : vertices(g))
    put(external_vertex_index_map, v, --vi);

  test_vertex_index_map_uniqueness(g, CGAL::parameters::vertex_index_map(external_vertex_index_map));

  // Read-only pmap
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor      halfedge_descriptor;
  typedef std::map<halfedge_descriptor, int>                            HalfedgeIndexMap;
  typedef CGAL::Non_mutable_property_map<halfedge_descriptor, int,
                                         HalfedgeIndexMap>              HalfedgeIdPropertyMap;

  int hi = 0;
  HalfedgeIndexMap him;
  HalfedgeIdPropertyMap external_halfedge_index_map(him);

  // this should complain that the map is not writable (commented because it does assert)
// CGAL::BGL::internal::initialize_halfedge_index_map(external_halfedge_index_map, g);

  // forced to initialize the underlying map
  for(halfedge_descriptor h : halfedges(g))
    him[h] = hi++;

  test_halfedge_index_map_uniqueness(g, CGAL::parameters::halfedge_index_map(external_halfedge_index_map));

  // Writable pmap
  typedef typename boost::graph_traits<Graph>::edge_descriptor          edge_descriptor;
  typedef boost::unordered_map<edge_descriptor, int>                    EdgeIndexMap;
  typedef CGAL::RW_property_map<edge_descriptor, int, EdgeIndexMap>     EdgeIdPropertyMap;

  EdgeIndexMap eim;
  EdgeIdPropertyMap external_edge_index_map(eim);
  CGAL::BGL::internal::initialize_edge_index_map(external_edge_index_map, g);

  test_edge_index_map_uniqueness(g, CGAL::parameters::edge_index_map(external_edge_index_map));

  // Just so face_index_map don't feel excluded
  typedef typename boost::graph_traits<Graph>::face_descriptor          face_descriptor;
  typedef std::map<face_descriptor, int>                                FaceIndexMap;
  typedef boost::const_associative_property_map<FaceIndexMap>           FaceIdPropertyMap;

  FaceIndexMap fim;
  FaceIdPropertyMap external_face_index_map(fim);

  // 'const_associative_pmap' has category 'lvalue_property_map_tag' but it's not writable
  // so below should complain (commented because it does assert)
//  CGAL::BGL::internal::initialize_face_index_map(external_face_index_map, g);

  // gotta initialize the underlying map
  int fi = 0;
  for(face_descriptor f : faces(g))
    fim[f] = fi++;

  test_face_index_map_uniqueness(g, CGAL::parameters::face_index_map(external_face_index_map));
}

template <typename Graph>
void test_all_index_maps_const(const Graph& g)
{
#ifdef CGAL_TEST_PROPERTIES_DEBUG
  std::cout << " ----------------------------  Const graph tests" << std::endl;
#endif

  test_internal_index_maps_const(g);
  test_initialized_index_maps_const(g);
}

///////////////////////////////////// non-const ////////////////////////////////////////////////////

template <typename Graph>
void test_internal_index_maps(Graph& g)
{
  test_uniqueness(g, vertices(g),  get(boost::vertex_index, g));
  test_uniqueness(g, halfedges(g), get(boost::halfedge_index, g));
  test_uniqueness(g, edges(g) ,    get(boost::edge_index, g));
  test_uniqueness(g, faces(g),     get(boost::face_index, g));
}

template <typename Graph>
void test_initialized_index_maps(Graph& g)
{
  typedef typename CGAL::GetInitializedVertexIndexMap<Graph>::type      VIM;
  VIM ivim = CGAL::get_initialized_vertex_index_map(g);
  test_uniqueness(g, vertices(g), ivim);

  typedef typename CGAL::GetInitializedHalfedgeIndexMap<Graph>::type    HIM;
  HIM ihim = CGAL::get_initialized_halfedge_index_map(g);
  test_uniqueness(g, halfedges(g), ihim);

  typedef typename CGAL::GetInitializedEdgeIndexMap<Graph>::type        EIM;
  EIM ieim = CGAL::get_initialized_edge_index_map(g);
  test_uniqueness(g, edges(g), ieim);

  typedef typename CGAL::GetInitializedFaceIndexMap<Graph>::type        FIM;
  FIM ifim = CGAL::get_initialized_face_index_map(g);
  test_uniqueness(g, faces(g), ifim);

  // Passing an index map via NP
  typedef typename boost::graph_traits<Graph>::vertex_descriptor        vertex_descriptor;
  typedef std::map<vertex_descriptor, int>                              VertexIndexMap;
  typedef boost::associative_property_map<VertexIndexMap>               VertexIdPropertyMap; // lvalue_pmap

  int vi = static_cast<int>(num_vertices(g));
  VertexIndexMap vim;
  VertexIdPropertyMap external_vertex_index_map(vim);
  for(vertex_descriptor v : vertices(g))
    put(external_vertex_index_map, v, --vi);

  test_vertex_index_map_uniqueness(g, CGAL::parameters::vertex_index_map(external_vertex_index_map));

  // Read-only pmap
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor      halfedge_descriptor;
  typedef std::map<halfedge_descriptor, int>                            HalfedgeIndexMap;
  typedef CGAL::Non_mutable_property_map<halfedge_descriptor, int,
                                         HalfedgeIndexMap>              HalfedgeIdPropertyMap;

  int hi = 0;
  HalfedgeIndexMap him;
  HalfedgeIdPropertyMap external_halfedge_index_map(him);

  // this should complain that the map is not writable (commented because it does assert)
// CGAL::BGL::internal::initialize_halfedge_index_map(external_halfedge_index_map, g);

  // forced to initialize the underlying map
  for(halfedge_descriptor h : halfedges(g))
    him[h] = hi++;

  test_halfedge_index_map_uniqueness(g, CGAL::parameters::halfedge_index_map(external_halfedge_index_map));

  // Writable pmap
  typedef typename boost::graph_traits<Graph>::edge_descriptor          edge_descriptor;
  typedef boost::unordered_map<edge_descriptor, int>                    EdgeIndexMap;
  typedef CGAL::RW_property_map<edge_descriptor, int, EdgeIndexMap>     EdgeIdPropertyMap;

  EdgeIndexMap eim;
  EdgeIdPropertyMap external_edge_index_map(eim);
  CGAL::BGL::internal::initialize_edge_index_map(external_edge_index_map, g);

  test_edge_index_map_uniqueness(g, CGAL::parameters::edge_index_map(external_edge_index_map));

  // Just so face_index_map don't feel excluded
  typedef typename boost::graph_traits<Graph>::face_descriptor          face_descriptor;
  typedef std::map<face_descriptor, int>                                FaceIndexMap;
  typedef boost::const_associative_property_map<FaceIndexMap>           FaceIdPropertyMap;

  FaceIndexMap fim;
  FaceIdPropertyMap external_face_index_map(fim);

  // 'const_associative_pmap' has category 'lvalue_property_map_tag' but it's not writable
  // so below should complain (commented because it does assert)
//  CGAL::BGL::internal::initialize_face_index_map(external_face_index_map, g);

  // gotta initialize the underlying map
  int fi = 0;
  for(face_descriptor f : faces(g))
    fim[f] = fi++;

  test_face_index_map_uniqueness(g, CGAL::parameters::face_index_map(external_face_index_map));
}

template <typename Graph>
void test_all_index_maps(Graph& g)
{
#ifdef CGAL_TEST_PROPERTIES_DEBUG
  std::cout << " ---------------------------- Non-const graph tests" << std::endl;
#endif

  test_internal_index_maps(g);
  test_initialized_index_maps(g);
}

template <typename Graph>
void test_graph(Graph& g)
{
#ifdef CGAL_TEST_PROPERTIES_DEBUG
  std::cout << "Graph has:" << std::endl
            << "\t" << num_vertices(g) << " vertices (actual: " << vertices(g).size() << ")" << std::endl
            << "\t" << num_halfedges(g) << " halfedges (actual: " << halfedges(g).size() << ")" << std::endl
            << "\t" << num_edges(g) << " edges (actual: " << edges(g).size() << ")" << std::endl
            << "\t" << num_faces(g) << " faces (actual: " << faces(g).size() << ")" << std::endl;
#endif

  test_all_index_maps(g);
  test_all_index_maps_const(g);
}

void test_poly(Polyhedron& g)
{
  test_graph(g);

  test_uniqueness(g, edges(g) ,    get(boost::edge_external_index, g));
  test_uniqueness(g, vertices(g),  get(boost::vertex_external_index, g));
  test_uniqueness(g, faces(g),     get(boost::face_external_index, g));
  test_uniqueness(g, halfedges(g), get(boost::halfedge_external_index, g));
}

int main(int, char**)
{
  std::cout << "testing Polyhedron\n";
  std::vector<Polyhedron> polys = poly_data();
  for(Polyhedron& p : polys)
    test_poly(p);

  std::cout << "testing Linear_cell_complex\n";
  std::vector<LCC> lccs = lcc_data();
  for(LCC& p : lccs)
    test_graph(p);

  std::cout << "testing Surface_mesh\n";
  std::vector<SM> sms = sm_data();
  for(SM& sm : sms)
  {
    assert(!CGAL::is_empty(sm));

    // Add some garbage
    CGAL::Euler::join_vertex(*(halfedges(sm).begin()), sm);

    test_graph(sm);

    // Test on a mesh with no internal index maps
    Seam_edge_pmap seam_edges = sm.add_property_map<SM::Edge_index, bool>("e:on_seam", false).first;
    Seam_vertex_pmap seam_vertices = sm.add_property_map<SM::Vertex_index, bool>("v:on_seam", false).first;
    Seam_mesh seam_mesh(sm, seam_edges, seam_vertices);

    test_initialized_index_maps(seam_mesh);
    test_initialized_index_maps_const(seam_mesh);
  }

#if defined(CGAL_USE_OPENMESH)
  std::cout << "testing OpenMesh\n";
  std::vector<OMesh> omeshs = omesh_data();
  for(OMesh& p : omeshs)
    test_graph(p);
#endif

  std::cout << "testing Triangulations\n";

  Triangulation_2 t2 = t2_data();
  test_graph(t2);

  Delaunay_triangulation_2 dt2 = dt2_data();
  test_graph(dt2);

  Regular_triangulation_2 rt2 = rt2_data();
  test_graph(rt2);

  Constrained_triangulation_2 ct2 = ct2_data();
  test_graph(ct2);

  Constrained_Delaunay_triangulation_2 cdt2 = cdt2_data();
  test_graph(cdt2);

  CDT_P2 cdtp2 = cdtp2_data();
  test_graph(cdtp2);

  Triangulation_hierarchy_2 t2h = t2h_data();
  test_graph(t2h);

  // no dynamic pmaps in triangulations (yet)
//  Triangulation_no_id_2 t2_no_id = t2_no_id_data();
//  test_initialized_index_maps(t2_no_id);
//  test_initialized_index_maps_const(t2_no_id);

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
