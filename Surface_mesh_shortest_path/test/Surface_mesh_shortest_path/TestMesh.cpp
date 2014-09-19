#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Surface_mesh_shortest_path/Surface_mesh_shortest_path_traits.h>
#include <CGAL/Surface_mesh_shortest_path/Surface_mesh_shortest_path.h>
#include <CGAL/Surface_mesh_shortest_path/function_objects.h>
#include <CGAL/Surface_mesh_shortest_path/barycentric.h>
#include <CGAL/Surface_mesh_shortest_path/internal/misc_functions.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

#include <boost/algorithm/string.hpp>

#include <CGAL/Random.h>

#include <CGAL/test_util.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <utility>
#include <cstdlib>
#include <cmath>

using namespace CGAL::test;

template <class Kernel>
struct TestMeshProgramInstance
{

//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Polyhedron_3> Traits;
typedef typename Traits::Barycentric_coordinate Barycentric_coordinate;
typedef typename Traits::FT FT;
typedef typename Traits::Point_3 Point_3;
typedef typename Traits::Point_2 Point_2;
typedef typename Traits::Triangle_3 Triangle_3;
typedef typename Traits::Triangle_2 Triangle_2;
typedef typename Traits::Segment_2 Segment_2;
typedef boost::graph_traits<Polyhedron_3> GraphTraits;
typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
typedef typename GraphTraits::vertex_iterator vertex_iterator;
typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
typedef typename GraphTraits::halfedge_iterator halfedge_iterator;
typedef typename GraphTraits::face_descriptor face_descriptor;
typedef typename GraphTraits::face_iterator face_iterator;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
typedef typename Surface_mesh_shortest_path::Face_location Face_location;
typedef typename boost::property_map<Polyhedron_3, CGAL::vertex_point_t>::type VPM;
typedef typename boost::property_map<typename Traits::FaceListGraph, boost::vertex_index_t>::type VIM;
typedef typename boost::property_map<typename Traits::FaceListGraph, boost::edge_index_t>::type EIM;
typedef typename boost::property_map<typename Traits::FaceListGraph, boost::halfedge_index_t>::type HIM;
typedef typename boost::property_map<typename Traits::FaceListGraph, boost::face_index_t>::type FIM;


TestMeshProgramInstance()
{
  randomSeed = 2681972;
  debugMode = false;
  randomizer = NULL;
  numVertices = 0;
}

size_t randomSeed;
size_t numIterations;
std::string meshName;
bool debugMode;

CGAL::Random* randomizer;
size_t numVertices;

Face_location next_location(Surface_mesh_shortest_path& shortestPath, Polyhedron_3& polyhedron, const std::vector<vertex_descriptor>& vertices)
{
  typename Traits::Construct_barycentric_coordinate construct_barycentric_coordinate;
  
  if (randomizer)
  {
    return shortestPath.face_location(vertices[randomizer->get_int(0, vertices.size())]);
  }
  else
  {
    std::string type;
      
    std::cin >> type;
    
    boost::algorithm::to_lower(type);
    
    if (type == "v")
    {
      size_t x;
      std::cin >> x;
      
      return shortestPath.face_location(vertices[x]);
    }
    else if (type == "e")
    {
      size_t x, y;
      double alpha;
      std::cin >> x >> y >> alpha;
      std::pair<halfedge_descriptor, bool> he = CGAL::halfedge(vertices[x], vertices[y], polyhedron);
      assert(he.second);
      return shortestPath.face_location(he.first, FT(alpha));
    }
    else if (type == "f")
    {
      size_t x, y;
      double alpha0, alpha1, alpha2;
      std::cin >> x >> y >> alpha0 >> alpha1 >> alpha2;
      std::pair<halfedge_descriptor, bool> he = CGAL::halfedge(vertices[x], vertices[y], polyhedron);
      return Face_location(CGAL::face(he.first, polyhedron), construct_barycentric_coordinate(FT(alpha0), FT(alpha1), FT(alpha2)));
    }
    
    return Face_location(GraphTraits::null_face(), construct_barycentric_coordinate(FT(0.0), FT(0.0), FT(0.0)));
  }
}

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
  Traits traits;

  Polyhedron_3 polyhedron;
  std::ifstream in(meshName.c_str());
  
  CGAL::set_halfedgeds_items_id(polyhedron);
  
  in >> polyhedron;
  
  in.close();
  
  numVertices = boost::num_vertices(polyhedron);
  
  VIM vertexIndexMap(get(boost::vertex_index, polyhedron));
  HIM halfedgeIndexMap(get(boost::halfedge_index, polyhedron));
  FIM faceIndexMap(get(boost::face_index, polyhedron));
  
  vertex_iterator verticesStart;
  vertex_iterator verticesEnd;
  
  std::vector<vertex_descriptor> vertices;
  
  boost::tie(verticesStart, verticesEnd) = boost::vertices(polyhedron);

  for (vertex_iterator it = verticesStart; it != verticesEnd; ++it)
  {
    vertices.push_back(*it);
  }
  
  face_iterator facesStart;
  face_iterator facesEnd;
  
  std::vector<face_descriptor> faces;
  
  boost::tie(facesStart, facesEnd) = CGAL::faces(polyhedron);
  
  for (face_iterator it = facesStart; it != facesEnd; ++it)
  {
    faces.push_back(*it);
  }

  Surface_mesh_shortest_path startToEndShortestPaths(polyhedron, traits);
  startToEndShortestPaths.m_debugOutput = debugMode;
  
  Surface_mesh_shortest_path endToStartShortestPaths(polyhedron, traits);
  endToStartShortestPaths.m_debugOutput = debugMode;

  std::cout << "Mesh: " << meshName << " " << boost::num_vertices(polyhedron) << " " << CGAL::num_faces(polyhedron) << " " << CGAL::num_halfedges(polyhedron) << std::endl;

  std::cout << std::setprecision(20);
  
  for (size_t i = 0; i < numIterations; ++i)
  {
    Face_location startLocation = next_location(startToEndShortestPaths, polyhedron, vertices);
    Face_location endLocation = next_location(endToStartShortestPaths, polyhedron, vertices);

    std::cout << "STE(location): " << faceIndexMap[startLocation.first] << " , " << startLocation.second[0] << " " << startLocation.second[1] << " " << startLocation.second[2] << " " << std::endl;
    std::cout << "ETS(location): " << faceIndexMap[endLocation.first] << " , " << endLocation.second[0] << " " << endLocation.second[1] << " " << endLocation.second[2] << " " << std::endl;

    std::cout << "Here" << __LINE__ << std::endl;
    
    startToEndShortestPaths.construct_sequence_tree(startLocation.first, startLocation.second);

    CGAL::test::Edge_sequence_collector<Traits> startToEndCollector(vertexIndexMap, halfedgeIndexMap, faceIndexMap);
    
    std::cout << "Here" << __LINE__ << std::endl;
    
    startToEndShortestPaths.shortest_path_sequence_to_source_points(endLocation.first, endLocation.second, startToEndCollector);
    
    std::cout << "Here" << __LINE__ << std::endl;

    FT startToEnd = startToEndShortestPaths.shortest_distance_to_source_points(endLocation.first, endLocation.second).first;
    
    std::cout << "Here" << __LINE__ << std::endl;

    endToStartShortestPaths.construct_sequence_tree(endLocation.first, endLocation.second);

    std::cout << "Here" << __LINE__ << std::endl;
    
    CGAL::test::Edge_sequence_collector<Traits> endToStartCollector(vertexIndexMap, halfedgeIndexMap, faceIndexMap);
    //endToStartShortestPaths.shortest_path_sequence_to_source_points(vertices[401], endToStartCollector);
    endToStartShortestPaths.shortest_path_sequence_to_source_points(startLocation.first, startLocation.second, endToStartCollector);
    
    //std::cout << "Weird: " << endToStartShortestPaths.shortest_distance_to_source_points(vertices[401]) << std::endl;

    FT endToStart = endToStartShortestPaths.shortest_distance_to_source_points(startLocation.first, startLocation.second).first;
    //
    
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
    Surface_mesh_shortest_path* pathStructures[2] = { &startToEndShortestPaths, &endToStartShortestPaths };
    CGAL::test::Edge_sequence_collector<Traits>* collectors[2] = { &startToEndCollector, &endToStartCollector };
    
    for (size_t d = 0; d < 2; ++d)
    {
      for (size_t j = 0; j < collectors[d]->m_sequence.size(); ++j)
      {
        std::cout << names[d] << "(sequence:" << j << "): ";
        
        CGAL::test::Sequence_item<Traits>& seqItem = collectors[d]->m_sequence[j];
        
        if (seqItem.type == CGAL::test::SEQUENCE_ITEM_EDGE)
        {
          std::cout << vertexIndexMap[CGAL::source(seqItem.halfedge, polyhedron)] << " , " << vertexIndexMap[CGAL::target(seqItem.halfedge, polyhedron)] << " : " << seqItem.edgeAlpha << std::endl;
        }
        else if (seqItem.type == CGAL::test::SEQUENCE_ITEM_VERTEX)
        {
          std::cout << vertexIndexMap[seqItem.vertex] << " , Distance: " << pathStructures[d]->shortest_distance_to_source_points(seqItem.vertex).first << std::endl;
        }
      }
    }
  }
}

};

int main(int argc, char** argv)
{
  TestMeshProgramInstance<CGAL::Exact_predicates_inexact_constructions_kernel> programInstance;

  if (argc > 4)
  {
    programInstance.randomSeed = std::atoi(argv[4]);
    programInstance.randomizer = new CGAL::Random(programInstance.randomSeed);
  }
  
  if (argc > 3)
  {
    programInstance.debugMode = (bool)std::atoi(argv[3]);
  }
  
  if (argc > 2)
  {
    programInstance.numIterations = std::atoi(argv[2]);
  }
  
  if (argc == 1)
  {
    std::cout << "Usage: " << argv[0] << " <mesh.off> <numIterations> <debugMode> <randseed>" << std::endl;
    return 0;
  }
  
  programInstance.meshName = std::string(argv[1]);
    
  programInstance.test_mesh_function();
  
  return 0;
}

// Hack to trick cgal_test_with_cmake into using this file even without a main
// int main(int argc, char** argv)
// {
//   return 0;
// }
