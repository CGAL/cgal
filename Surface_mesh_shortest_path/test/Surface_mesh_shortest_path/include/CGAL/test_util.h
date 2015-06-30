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
  SEQUENCE_ITEM_FACE
};

template <class Traits>
struct Sequence_item
{
  typedef typename Traits::Triangle_mesh Triangle_mesh;
  typedef typename Traits::FT FT;
  typedef typename Traits::Barycentric_coordinate Barycentric_coordinate;
  typedef typename boost::graph_traits<Triangle_mesh> Graph_traits;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  
  Sequence_item_type type;
  size_t index;
  Barycentric_coordinate faceAlpha;
  FT edgeAlpha;
  
  halfedge_descriptor halfedge;
  vertex_descriptor vertex;
  face_descriptor face;
};

template <class Traits, 
  class VIM = typename boost::property_map<typename Traits::Triangle_mesh, boost::vertex_index_t>::type,
  class HIM = typename boost::property_map<typename Traits::Triangle_mesh, boost::halfedge_index_t>::type,
  class FIM = typename boost::property_map<typename Traits::Triangle_mesh, boost::face_index_t>::type>
struct Edge_sequence_collector
{
  typedef typename Traits::Triangle_mesh Triangle_mesh;
  typedef typename Traits::FT FT;
  typedef typename Traits::Barycentric_coordinate Barycentric_coordinate;
  typedef VIM VertexIndexMap;
  typedef HIM HalfedgeIndexMap;
  typedef FIM FaceIndexMap;
  typedef typename boost::graph_traits<Triangle_mesh> Graph_traits;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename Graph_traits::face_descriptor face_descriptor;

  VertexIndexMap m_vertexIndexMap;
  HalfedgeIndexMap m_halfedgeIndexMap;
  FaceIndexMap m_faceIndexMap;
  
  std::vector<Sequence_item<Traits> > m_sequence;
  
  Edge_sequence_collector(Triangle_mesh& g)
    : m_vertexIndexMap(get(boost::vertex_index, g))
    , m_halfedgeIndexMap(get(CGAL::halfedge_index, g))
    , m_faceIndexMap(get(CGAL::face_index, g))
  {
  }

  Edge_sequence_collector(VertexIndexMap& vertexIndexMap, HalfedgeIndexMap& halfedgeIndexMap, FaceIndexMap& faceIndexMap)
    : m_vertexIndexMap(vertexIndexMap)
    , m_halfedgeIndexMap(halfedgeIndexMap)
    , m_faceIndexMap(faceIndexMap)
  {
  }
  
  void operator()(halfedge_descriptor he, FT alpha)
  {
    Sequence_item<Traits> item;
    item.type = SEQUENCE_ITEM_EDGE;
    item.index = m_halfedgeIndexMap[he];
    item.edgeAlpha = alpha;
    item.halfedge = he;
    m_sequence.push_back(item);
  }
  
  void operator()(vertex_descriptor v)
  {
    Sequence_item<Traits> item;
    item.type = SEQUENCE_ITEM_VERTEX;
    item.index = m_vertexIndexMap[v];
    item.vertex = v;
    m_sequence.push_back(item);
  }
  
  void operator()(face_descriptor f, Barycentric_coordinate alpha)
  {
    Sequence_item<Traits> item;
    item.type = SEQUENCE_ITEM_FACE;
    item.index = m_faceIndexMap[f];
    item.faceAlpha = alpha;
    item.face = f;
    m_sequence.push_back(item);
  }
};

template <class Traits>
typename Traits::Barycentric_coordinate random_coordinate(CGAL::Random& rand)
{
  typedef typename Traits::FT FT;
  typename Traits::Construct_barycentric_coordinate construct_barycentric_coordinate;
  FT u = rand.uniform_real(FT(0.0), FT(1.0));
  FT v = rand.uniform_real(FT(0.0), FT(FT(1.0) - u));
  return construct_barycentric_coordinate(u, v, FT(FT(1.0) - u - v));
}

template <class FT>
FT squared(FT in)
{
  return in * in;
}

template <class Triangle_mesh>
typename Triangle_mesh::Halfedge_handle make_regular_tetrahedron(Triangle_mesh& out)
{
  typedef typename Triangle_mesh::Traits::FT FT;
  
  FT rsqrt2 = FT(1.0) / CGAL::sqrt(FT(2.0));
  out.clear();
  typename Triangle_mesh::Halfedge_handle result = out.make_tetrahedron(
  typename Triangle_mesh::Point_3(FT(1.0), FT(0.0), -rsqrt2),
  typename Triangle_mesh::Point_3(-FT(1.0), FT(0.0), -rsqrt2),
  typename Triangle_mesh::Point_3(FT(0.0), FT(1.0), rsqrt2),
  typename Triangle_mesh::Point_3(FT(0.0), -FT(1.0), rsqrt2));

  return result;
}

template <class Triangle_mesh>
size_t face_vertex_index(typename boost::graph_traits<Triangle_mesh>::face_descriptor face, typename boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex, Triangle_mesh& g)
{
  size_t index = 0;
  
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor halfedge_descriptor;
  
  halfedge_descriptor currentEdge(CGAL::halfedge(face, g));
  halfedge_descriptor startEdge = currentEdge;
  
  do
  {
    if (CGAL::source(currentEdge, g) == vertex)
    {
      return index;
    }
    
    ++index;
    currentEdge = CGAL::next(currentEdge, g);
  }
  while (currentEdge != startEdge);
  
  return index;
}

} // namespace util

} // namespace CGAL

#endif // CGAL_TEST_UTIL_H
