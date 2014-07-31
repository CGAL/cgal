#ifndef CGAL_TEST_UTIL_H
#define CGAL_TEST_UTIL_H

#include <algorithm>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Random.h>

namespace CGAL {

namespace test {

enum Sequence_item_type
{
  SEQUENCE_ITEM_VERTEX,
  SEQUENCE_ITEM_EDGE,
  SEQUENCE_ITEM_FACE,
};

template <class Traits>
struct Sequence_item
{
  typedef typename Traits::FaceGraph FaceGraph;
  typedef typename Traits::FT FT;
  typedef typename Traits::Barycentric_coordinate Barycentric_coordinate;
  typedef typename boost::graph_traits<FaceGraph> GraphTraits;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  
  Sequence_item_type type;
  size_t index;
  Barycentric_coordinate faceAlpha;
  FT edgeAlpha;
  
  halfedge_descriptor halfedge;
  vertex_descriptor vertex;
  face_descriptor face;
};

template <class Traits, 
  class VIM = typename boost::property_map<typename Traits::FaceGraph, CGAL::vertex_external_index_t>::type,
  class HIM = typename boost::property_map<typename Traits::FaceGraph, CGAL::halfedge_external_index_t>::type,
  class FIM = typename boost::property_map<typename Traits::FaceGraph, CGAL::face_external_index_t>::type>
struct Edge_sequence_collector
{
  typedef typename Traits::FaceGraph FaceGraph;
  typedef typename Traits::FT FT;
  typedef typename Traits::Barycentric_coordinate Barycentric_coordinate;
  typedef VIM VertexIndexMap;
  typedef HIM HalfedgeIndexMap;
  typedef FIM FaceIndexMap;
  typedef typename boost::graph_traits<FaceGraph> GraphTraits;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;

  VertexIndexMap m_vertexIndexMap;
  HalfedgeIndexMap m_halfedgeIndexMap;
  FaceIndexMap m_faceIndexMap;
  
  std::vector<Sequence_item<Traits> > m_sequence;
  
  Edge_sequence_collector(FaceGraph& p)
    : m_vertexIndexMap(CGAL::get(boost::vertex_external_index, p))
    , m_halfedgeIndexMap(CGAL::get(CGAL::halfedge_external_index, p))
    , m_faceIndexMap(CGAL::get(CGAL::face_external_index, p))
  {
  }

  Edge_sequence_collector(VertexIndexMap& vertexIndexMap, HalfedgeIndexMap& halfedgeIndexMap, FaceIndexMap& faceIndexMap)
    : m_vertexIndexMap(vertexIndexMap)
    , m_halfedgeIndexMap(halfedgeIndexMap)
    , m_faceIndexMap(faceIndexMap)
  {
  }
  
  void edge(halfedge_descriptor he, FT alpha)
  {
    Sequence_item<Traits> item;
    item.type = SEQUENCE_ITEM_EDGE;
    item.index = m_halfedgeIndexMap[he];
    item.edgeAlpha = alpha;
    item.halfedge = he;
    m_sequence.push_back(item);
  }
  
  void vertex(vertex_descriptor v)
  {
    Sequence_item<Traits> item;
    item.type = SEQUENCE_ITEM_VERTEX;
    item.index = m_vertexIndexMap[v];
    item.vertex = v;
    m_sequence.push_back(item);
  }
  
  void face(face_descriptor f, Barycentric_coordinate alpha)
  {
    Sequence_item<Traits> item;
    item.type = SEQUENCE_ITEM_FACE;
    item.index = m_faceIndexMap[f];
    item.faceAlpha = alpha;
    item.face = f;
    m_sequence.push_back(item);
  }
};

template <class FT, class B>
B random_coordinate(CGAL::Random& rand)
{
  FT u = rand.uniform_01<FT>();
  FT v = rand.uniform_real(FT(0.0), FT(1.0) - u);
  return B(u, v, FT(1.0) - u - v);
}

template <class FT>
FT squared(FT in)
{
  return in * in;
}

template<class FaceGraph>
struct Plane_from_facet {
  typedef typename FaceGraph::Plane_3 Plane_3;
  typedef typename FaceGraph::Facet Facet;
  typedef typename FaceGraph::Halfedge_handle Halfedge_handle;

  Plane_3 operator()(Facet& f) {
      Halfedge_handle h = f.halfedge();
      return Plane_3( h->vertex()->point(),
                                    h->next()->vertex()->point(),
                                    h->opposite()->vertex()->point());
  }
};

template <class FaceGraph>
void construct_polyhedron_planes(FaceGraph& out)
{
  std::transform( out.facets_begin(), out.facets_end(), out.planes_begin(), Plane_from_facet<FaceGraph>());
}

template <class FaceGraph>
typename FaceGraph::Halfedge_handle make_regular_tetrahedron(FaceGraph& out)
{
  typedef typename FaceGraph::Traits::FT FT;
  
  FT rsqrt2 = FT(1.0) / CGAL::sqrt(FT(2.0));
  out.clear();
  typename FaceGraph::Halfedge_handle result = out.make_tetrahedron(
    typename FaceGraph::Point_3(FT(1.0), FT(0.0), -rsqrt2),
    typename FaceGraph::Point_3(-FT(1.0), FT(0.0), -rsqrt2),
    typename FaceGraph::Point_3(FT(0.0), FT(1.0), rsqrt2),
    typename FaceGraph::Point_3(FT(0.0), -FT(1.0), rsqrt2));
  construct_polyhedron_planes(out);
  return result;
}

template <class FaceGraph>
size_t face_vertex_index(typename boost::graph_traits<FaceGraph>::face_descriptor face, typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex, FaceGraph& P)
{
  size_t index = 0;
  
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  
  halfedge_descriptor currentEdge(CGAL::halfedge(face, P));
  halfedge_descriptor startEdge = currentEdge;
  
  do
  {
    if (CGAL::source(currentEdge, P) == vertex)
    {
      return index;
    }
    
    ++index;
    currentEdge = CGAL::next(currentEdge, P);
  }
  while (currentEdge != startEdge);
  
  return index;
}

} // namespace util

} // namespace CGAL

#endif // CGAL_TEST_UTIL_H