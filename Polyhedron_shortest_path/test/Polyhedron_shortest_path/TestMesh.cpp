// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path_traits.h>
#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path.h>
#include <CGAL/Polyhedron_shortest_path/Internal/function_objects.h>
#include <CGAL/Polyhedron_shortest_path/Internal/Barycentric.h>
#include <CGAL/Polyhedron_shortest_path/Internal/misc_functions.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

#include <CGAL/Random.h>

#include <CGAL/test_util.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <utility>
#include <cstdlib>
#include <cmath>

/*
#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;
*/

enum Sequence_item_type
{
  SEQUENCE_ITEM_VERTEX,
  SEQUENCE_ITEM_EDGE,
  SEQUENCE_ITEM_FACE,
};

template <class Traits>
struct Sequence_item
{
  typedef typename Traits::Polyhedron Polyhedron;
  typedef typename Traits::FT FT;
  typedef typename Traits::Barycentric_coordinate Barycentric_coordinate;
  typedef typename boost::graph_traits<Polyhedron> GraphTraits;
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
  class VIM = typename boost::property_map<typename Traits::Polyhedron, CGAL::vertex_external_index_t>::type,
  class HIM = typename boost::property_map<typename Traits::Polyhedron, CGAL::halfedge_external_index_t>::type,
  class FIM = typename boost::property_map<typename Traits::Polyhedron, CGAL::face_external_index_t>::type>
struct Edge_sequence_collector
{
  typedef typename Traits::Polyhedron Polyhedron;
  typedef typename Traits::FT FT;
  typedef typename Traits::Barycentric_coordinate Barycentric_coordinate;
  typedef VIM VertexIndexMap;
  typedef HIM HalfedgeIndexMap;
  typedef FIM FaceIndexMap;
  typedef typename boost::graph_traits<Polyhedron> GraphTraits;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;

  VertexIndexMap m_vertexIndexMap;
  HalfedgeIndexMap m_halfedgeIndexMap;
  FaceIndexMap m_faceIndexMap;
  
  std::vector<Sequence_item<Traits> > m_sequence;
  
  Edge_sequence_collector(Polyhedron& p)
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

size_t randomSeed = 2681972;
size_t numIterations;
std::string meshName;
bool debugMode = false;

CGAL::Random* randomizer = NULL;
size_t numVertices = 0;

size_t next_vertex()
{
  if (randomizer)
  {
    return randomizer->get_int(0, numVertices);
  }
  else
  {
    size_t x;
    std::cin >> x;
    return x;
  }
}

void test_mesh_function()
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;
  typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Polyhedron_3> Traits;
  typedef Traits::Barycentric_coordinate Barycentric_coordinate;
  typedef Traits::FT FT;
  typedef Traits::Point_3 Point_3;
  typedef Traits::Point_2 Point_2;
  typedef Traits::Triangle_3 Triangle_3;
  typedef Traits::Triangle_2 Triangle_2;
  typedef Traits::Segment_2 Segment_2;
  typedef boost::graph_traits<Polyhedron_3> GraphTraits;
  typedef GraphTraits::vertex_descriptor vertex_descriptor;
  typedef GraphTraits::vertex_iterator vertex_iterator;
  typedef GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef GraphTraits::halfedge_iterator halfedge_iterator;
  typedef GraphTraits::face_descriptor face_descriptor;
  typedef GraphTraits::face_iterator face_iterator;
  typedef CGAL::Polyhedron_shortest_path<Traits> Polyhedron_shortest_path;
  typedef boost::property_map<Polyhedron_3, CGAL::vertex_point_t>::type VPM;
  typedef boost::property_map<typename Traits::Polyhedron, boost::vertex_external_index_t>::type VIM;
  typedef boost::property_map<typename Traits::Polyhedron, boost::edge_external_index_t>::type EIM;
  typedef boost::property_map<typename Traits::Polyhedron, CGAL::halfedge_external_index_t>::type HIM;
  typedef boost::property_map<typename Traits::Polyhedron, CGAL::face_external_index_t>::type FIM;
  
  Traits traits;
  
  Polyhedron_3 P;
  std::ifstream in(meshName.c_str());
  
  in >> P;
  
  in.close();
  
  numVertices = boost::num_vertices(P);
  
  VIM vertexIndexMap(CGAL::get(boost::vertex_external_index, P));
  HIM halfedgeIndexMap(CGAL::get(CGAL::halfedge_external_index, P));
  FIM faceIndexMap(CGAL::get(CGAL::face_external_index, P));
  
  vertex_iterator verticesStart;
  vertex_iterator verticesEnd;
  
  std::vector<vertex_descriptor> vertices;
  
  boost::tie(verticesStart, verticesEnd) = boost::vertices(P);
  
  size_t currentId = 0;
  
  for (vertex_iterator it = verticesStart; it != verticesEnd; ++it)
  {
    vertices.push_back(*it);
    //it->id() = currentId;
    //++currentId;
  }
  
  face_iterator facesStart;
  face_iterator facesEnd;
  
  std::vector<face_descriptor> faces;
  
  boost::tie(facesStart, facesEnd) = CGAL::faces(P);
  
  for (face_iterator it = facesStart; it != facesEnd; ++it)
  {
    faces.push_back(*it);
  }

  Polyhedron_shortest_path startToEndShortestPaths(P, traits);
  startToEndShortestPaths.m_debugOutput = debugMode;
  
  Polyhedron_shortest_path endToStartShortestPaths(P, traits);
  endToStartShortestPaths.m_debugOutput = debugMode;

  std::cout << "Mesh: " << meshName << " " << boost::num_vertices(P) << " " << CGAL::num_faces(P) << " " << CGAL::num_halfedges(P) << std::endl;

  for (size_t i = 0; i < numIterations; ++i)
  {
    bool found = false;
    size_t startVertexIndex;
    size_t endVertexIndex;
    vertex_descriptor startVertex;
    vertex_descriptor endVertex;
    
    startVertexIndex = next_vertex();
    endVertexIndex = next_vertex();
    
    std::cout << "STE(index): " << startVertexIndex << std::endl;
    std::cout << "ETS(index): " << endVertexIndex << std::endl;
    
    startVertex = vertices[startVertexIndex];
    endVertex = vertices[endVertexIndex];
    
    startToEndShortestPaths.compute_shortest_paths(startVertex);

    Edge_sequence_collector<Traits> startToEndCollector(vertexIndexMap, halfedgeIndexMap, faceIndexMap);
    startToEndShortestPaths.shortest_path_sequence(endVertex, startToEndCollector);

    FT startToEnd = startToEndShortestPaths.shortest_distance_to_vertex(endVertex);

    endToStartShortestPaths.compute_shortest_paths(endVertex);

    Edge_sequence_collector<Traits> endToStartCollector(vertexIndexMap, halfedgeIndexMap, faceIndexMap);
    endToStartShortestPaths.shortest_path_sequence(startVertex, endToStartCollector);

    FT endToStart = endToStartShortestPaths.shortest_distance_to_vertex(startVertex);

    std::cout << "STE(distance): " << startToEnd << std::endl;
    std::cout << "ETS(distance): " << endToStart << std::endl;
    if (CGAL::abs(startToEnd - endToStart) > FT(0.000001))
    {
      std::cout << "STE/ETS: distanceerror!" << std::endl;
    }
    std::cout << "STE(size): " << startToEndCollector.m_sequence.size() << std::endl;
    std::cout << "ETS(size): " << endToStartCollector.m_sequence.size() << std::endl;
    if (startToEndCollector.m_sequence.size() != endToStartCollector.m_sequence.size())
    {
      std::cout << "STE/ETS: sizeerror!" << std::endl;
    }
    
    
    std::string names[2] = { "STE", "ETS" };
    Polyhedron_shortest_path* pathStructures[2] = { &startToEndShortestPaths, &endToStartShortestPaths };
    Edge_sequence_collector<Traits>* collectors[2] = { &startToEndCollector, &endToStartCollector };
    
    for (size_t d = 0; d < 2; ++d)
    {
      for (size_t j = 0; j < collectors[d]->m_sequence.size(); ++j)
      {
        std::cout << names[d] << "(sequence:" << j << "): ";
        
        Sequence_item<Traits>& seqItem = collectors[d]->m_sequence[j];
        
        if (seqItem.type == SEQUENCE_ITEM_EDGE)
        {
          std::cout << vertexIndexMap[CGAL::source(seqItem.halfedge, P)] << " , " << vertexIndexMap[CGAL::target(seqItem.halfedge, P)] << " : " << seqItem.edgeAlpha << std::endl;
        }
        else if (seqItem.type == SEQUENCE_ITEM_VERTEX)
        {
          std::cout << vertexIndexMap[seqItem.vertex] << " , Distance: " << pathStructures[d]->shortest_distance_to_vertex(seqItem.vertex) << std::endl;
        }
      }
    }
  }
}

int main(int argc, char** argv)
{
  if (argc > 4)
  {
    randomSeed = std::atoi(argv[4]);
    randomizer = new CGAL::Random(randomSeed);
  }
  
  if (argc > 3)
  {
    debugMode = (bool)std::atoi(argv[3]);
  }
  
  if (argc > 2)
  {
    numIterations = std::atoi(argv[2]);
  }
  
  if (argc == 1)
  {
    std::cout << "Usage: " << argv[0] << " <mesh.off> <numIterations> <debugMode> <randseed>" << std::endl;
    return 0;
  }
  
  meshName = std::string(argv[1]);
    
  test_mesh_function();
  
  return 0;
}

// Hack to trick cgal_test_with_cmake into using this file even without a main
// int main(int argc, char** argv)
// {
//   return 0;
// }