#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Random.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>


#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/Polyhedron_shortest_path.h>

enum Sequence_item_type
{
  SEQUENCE_ITEM_VERTEX,
  SEQUENCE_ITEM_EDGE,
  SEQUENCE_ITEM_FACE,
};

// Stores a single item in a shortest paths sequence (edge, vertex, crossed, or a face end-point).
template <class Traits>
struct Sequence_item
{
  typedef typename Traits::FaceListGraph FaceListGraph;
  typedef typename Traits::FT FT;
  typedef typename Traits::Barycentric_coordinate Barycentric_coordinate;
  typedef typename boost::graph_traits<FaceListGraph> GraphTraits;
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

// A model of PolyhedronShortestPathVisitor
template <class Traits, 
  class VIM = typename boost::property_map<typename Traits::FaceListGraph, boost::vertex_index_t>::type,
  class HIM = typename boost::property_map<typename Traits::FaceListGraph, boost::halfedge_index_t>::type,
  class FIM = typename boost::property_map<typename Traits::FaceListGraph, boost::face_index_t>::type>
struct Sequence_collector
{
  typedef typename Traits::FaceListGraph FaceListGraph;
  typedef typename Traits::FT FT;
  typedef typename Traits::Barycentric_coordinate Barycentric_coordinate;
  typedef VIM VertexIndexMap;
  typedef HIM HalfedgeIndexMap;
  typedef FIM FaceIndexMap;
  typedef typename boost::graph_traits<FaceListGraph> GraphTraits;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;

  VertexIndexMap m_vertexIndexMap;
  HalfedgeIndexMap m_halfedgeIndexMap;
  FaceIndexMap m_faceIndexMap;
  
  std::vector<Sequence_item<Traits> > m_sequence;
  
  Sequence_collector(FaceListGraph& g)
    : m_vertexIndexMap(CGAL::get(boost::vertex_index, g))
    , m_halfedgeIndexMap(CGAL::get(CGAL::halfedge_index, g))
    , m_faceIndexMap(CGAL::get(CGAL::face_index, g))
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

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;
typedef CGAL::Polyhedron_shortest_path_traits<Kernel, Polyhedron_3> Traits;
typedef CGAL::Polyhedron_shortest_path<Traits> Polyhedron_shortest_path;
typedef boost::graph_traits<Polyhedron_3> GraphTraits;
typedef GraphTraits::vertex_iterator vertex_iterator;
typedef GraphTraits::face_descriptor face_descriptor;
typedef GraphTraits::face_iterator face_iterator;

Traits::Barycentric_coordinate random_coordinate(CGAL::Random& rand)
{
  typename Traits::Construct_barycentric_coordinate construct_barycentric_coordinate;
  Traits::FT u = rand.uniform_real(Traits::FT(0.0), Traits::FT(1.0));
  Traits::FT v = rand.uniform_real(Traits::FT(0.0), Traits::FT(Traits::FT(1.0) - u));
  return construct_barycentric_coordinate(u, v, Traits::FT(Traits::FT(1.0) - u - v));
}

int main(int argc, char** argv)
{
  Polyhedron_3 polyhedron;
  
  std::ifstream inStream(argv[1]);
  
  inStream >> polyhedron;
  
  inStream.close();
  
  const size_t randSeed = argc > 2 ? std::atoi(argv[2]) : 8031760;
  CGAL::Random rand(randSeed);

  face_iterator facesStart, facesEnd;
  boost::tie(facesStart, facesEnd) = faces(polyhedron);
  
  std::vector<face_descriptor> faceList;
  
  for (face_iterator facesCurrent = facesStart; facesCurrent != facesEnd; ++facesCurrent)
  {
    faceList.push_back(*facesCurrent);
  }

  std::vector<std::pair<Polyhedron_shortest_path::Face_location, Polyhedron_shortest_path::Face_location > > locationPairs;

  Traits traits;
  Polyhedron_shortest_path shortestPaths(polyhedron, traits);
  
  Polyhedron_shortest_path::Face_location startLocation(faceList[rand.get_int(0, CGAL::num_faces(polyhedron))], random_coordinate(rand));
  Polyhedron_shortest_path::Face_location endLocation(faceList[rand.get_int(0, CGAL::num_faces(polyhedron))], random_coordinate(rand));
  
  shortestPaths.construct_sequence_tree(startLocation.first, startLocation.second);

  Sequence_collector<Traits> sequenceCollector(polyhedron);
  
  shortestPaths.shortest_path_sequence_to_source_points(endLocation.first, endLocation.second, sequenceCollector);
  
  for (size_t i = 0; i < sequenceCollector.m_sequence.size(); ++i)
  {
    Sequence_item<Traits>& item = sequenceCollector.m_sequence[i];
    
    switch (item.type)
    {
      case SEQUENCE_ITEM_VERTEX:
        std::cout << "#" << i << " : Vertex : " << item.vertex->id() << std::endl;
        break;
      case SEQUENCE_ITEM_EDGE:
        std::cout << "#" << i << " : Edge : " << item.halfedge->id() << " , (" << Traits::FT(1.0) - item.edgeAlpha << " , " << item.edgeAlpha << ")" << std::endl;
        break;
      case SEQUENCE_ITEM_FACE:
        std::cout << "#" << i << " : Face : " << item.face->id() << " , (" << item.faceAlpha[0] << " , " << item.faceAlpha[1] << " , " << item.faceAlpha[2] << ")" << std::endl;
        break;
    }
  }
  
  return 0;
}
