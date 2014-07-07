// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#include <CGAL/Interval_nt.h>
#include <iomanip>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path_traits.h>
#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path.h>
#include <CGAL/Polyhedron_shortest_path/Internal/function_objects.h>
#include <CGAL/Polyhedron_shortest_path/Internal/Barycentric.h>
#include <CGAL/Polyhedron_shortest_path/Internal/misc_functions.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/Random.h>

#include <CGAL/test_util.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <cmath>

#define BOOST_TEST_MODULE polyhedron_shortest_path_test_2
#include <boost/test/included/unit_test.hpp>

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

BOOST_AUTO_TEST_CASE( test_a_to_b_vs_b_t_a_distances )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  //typedef CGAL::Simple_cartesian<Auto_interval_nt> Kernel;
  typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
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
  
  CGAL::Random rand(2681972);
  
  std::string meshes[3] = { "data/elephant.off", "data/saddle_vertex_mesh.off", "data/anchor.off" };

  for (size_t meshId = 0; meshId < 3; ++meshId)
  {
    Polyhedron_3 P;
    std::ifstream in(meshes[meshId].c_str());
    
    in >> P;
    
    in.close();
    
    VIM vertexIndexMap(CGAL::get(boost::vertex_external_index, P));
    HIM halfedgeIndexMap(CGAL::get(CGAL::halfedge_external_index, P));
    FIM faceIndexMap(CGAL::get(CGAL::face_external_index, P));
    
    vertex_iterator verticesStart;
    vertex_iterator verticesEnd;
    
    std::vector<vertex_descriptor> vertices;
    
    boost::tie(verticesStart, verticesEnd) = boost::vertices(P);
    
    for (vertex_iterator it = verticesStart; it != verticesEnd; ++it)
    {
      vertices.push_back(*it);
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
    Polyhedron_shortest_path endToStartShortestPaths(P, traits);
    
    const size_t numTests = 15;
    
    std::cout << "Mesh: " << meshes[meshId] << std::endl;
    
    
    for (size_t i = 0; i < numTests; ++i)
    {
      size_t startVertexIndex = rand.get_int(0, vertices.size());
      size_t endVertexIndex = rand.get_int(0, vertices.size());
      
      vertex_descriptor startVertex = vertices[startVertexIndex];
      vertex_descriptor endVertex = vertices[endVertexIndex];

      //startToEndShortestPaths.m_debugOutput = true;
      
      startToEndShortestPaths.compute_shortest_paths(startVertex);

      Edge_sequence_collector<Traits> startToEndCollector(vertexIndexMap, halfedgeIndexMap, faceIndexMap);
      startToEndShortestPaths.shortest_path_sequence(endVertex, startToEndCollector);
      
      FT startToEnd = startToEndShortestPaths.shortest_distance_to_vertex(endVertex);
      
      //endToStartShortestPaths.m_debugOutput = true;
      
      endToStartShortestPaths.compute_shortest_paths(endVertex);
      
      Edge_sequence_collector<Traits> endToStartCollector(vertexIndexMap, halfedgeIndexMap, faceIndexMap);
      endToStartShortestPaths.shortest_path_sequence(startVertex, endToStartCollector);
      
      FT endToStart = endToStartShortestPaths.shortest_distance_to_vertex(startVertex);

      BOOST_CHECK_CLOSE(startToEnd, endToStart, FT(0.0000001));
      
      BOOST_CHECK_EQUAL(startToEndCollector.m_sequence.size(), endToStartCollector.m_sequence.size());
      
      if (startToEndCollector.m_sequence.size() == endToStartCollector.m_sequence.size())
      {
        if (startToEndCollector.m_sequence.size() > 3)
        {
          size_t k = startToEndCollector.m_sequence.size() - 2;
          
          for (size_t j = 0; j < endToStartCollector.m_sequence.size() - 1; ++j)
          {
            BOOST_CHECK_EQUAL(endToStartCollector.m_sequence[j].type, startToEndCollector.m_sequence[k].type);
            
            if (endToStartCollector.m_sequence[j].type == startToEndCollector.m_sequence[k].type)
            {
              switch (endToStartCollector.m_sequence[j].type)
              {
              case SEQUENCE_ITEM_VERTEX:
              case SEQUENCE_ITEM_FACE:
                BOOST_CHECK_EQUAL(endToStartCollector.m_sequence[j].index, startToEndCollector.m_sequence[k].index);
                break;
              case SEQUENCE_ITEM_EDGE:
                BOOST_CHECK_EQUAL(halfedgeIndexMap[endToStartCollector.m_sequence[j].halfedge], halfedgeIndexMap[CGAL::opposite(startToEndCollector.m_sequence[k].halfedge, P)]);
                break;
              }
            }
            else
            {
              Sequence_item<Traits>& steItem = startToEndCollector.m_sequence[k];
              Sequence_item<Traits>& etsItem = endToStartCollector.m_sequence[j];
              
              std::string names[2] = { "STE", "ETS" };
              Sequence_item<Traits> items[2] = { startToEndCollector.m_sequence[k], endToStartCollector.m_sequence[j] };
              Polyhedron_shortest_path* pathStructures[2] = { &startToEndShortestPaths, &endToStartShortestPaths };
              
              for (size_t d = 0; d < 2; ++d)
              {
                if (items[d].type == SEQUENCE_ITEM_EDGE)
                {
                  std::cout << "\t" << names[d] << "(edge): " << vertexIndexMap[CGAL::source(items[d].halfedge, P)] << " , " << vertexIndexMap[CGAL::target(items[d].halfedge, P)] << " : " << items[d].edgeAlpha << std::endl;
                }
                else if (items[d].type == SEQUENCE_ITEM_VERTEX)
                {
                  std::cout << "\t" << names[d] << "(vertex): " << vertexIndexMap[items[d].vertex] << " , Distance: " << pathStructures[d]->shortest_distance_to_vertex(items[d].vertex) << std::endl;
                }
              }
            }
            
            if (k > 0)
            {
              --k;
            }
          }
        }
      }
    }
    

    for (size_t i = 0; i < numTests; ++i)
    {
      size_t startFaceIndex = rand.get_int(0, faces.size());
      size_t endFaceIndex = rand.get_int(0, faces.size());
      
      face_descriptor startFace = faces[startFaceIndex];
      face_descriptor endFace = faces[endFaceIndex];
      
      Barycentric_coordinate startLocation = random_coordinate<FT, Barycentric_coordinate>(rand);
      Barycentric_coordinate endLocation = random_coordinate<FT, Barycentric_coordinate>(rand);

      //shortestPaths.m_debugOutput = true;
      
      startToEndShortestPaths.compute_shortest_paths(startFace, startLocation);

      //CGAL::Interval_nt<true> startToEnd = startToEndShortestPaths.shortest_distance_to_location_interval(endFace, endLocation);
      
      FT startToEnd = startToEndShortestPaths.shortest_distance_to_location(endFace, endLocation);
      
      Edge_sequence_collector<Traits> startToEndCollector(vertexIndexMap, halfedgeIndexMap, faceIndexMap);
      startToEndShortestPaths.shortest_path_sequence(endFace, endLocation, startToEndCollector);
      
      endToStartShortestPaths.compute_shortest_paths(endFace, endLocation);
      
      //CGAL::Interval_nt<true> endToStart = endToStartShortestPaths.shortest_distance_to_location_interval(startFace, startLocation);
      
      FT endToStart = endToStartShortestPaths.shortest_distance_to_location(startFace, startLocation);

      Edge_sequence_collector<Traits> endToStartCollector(vertexIndexMap, halfedgeIndexMap, faceIndexMap);
      endToStartShortestPaths.shortest_path_sequence(startFace, startLocation, endToStartCollector);
      
      //std::cout << std::setprecision(15) << std::endl;

      BOOST_CHECK_CLOSE(startToEnd, endToStart, FT(0.0000001));

      BOOST_CHECK_EQUAL(startToEndCollector.m_sequence.size(), endToStartCollector.m_sequence.size());
      
      if (startToEndCollector.m_sequence.size() > 3 && startToEndCollector.m_sequence.size() == endToStartCollector.m_sequence.size())
      {
      
        for (size_t j = 0; j < endToStartCollector.m_sequence.size() - 1; ++j)
        {
          size_t k = endToStartCollector.m_sequence.size() - j - 2;
          
          BOOST_CHECK_EQUAL(endToStartCollector.m_sequence[j].type, startToEndCollector.m_sequence[k].type);
          
          if (endToStartCollector.m_sequence[j].type == startToEndCollector.m_sequence[k].type)
          {
            switch (endToStartCollector.m_sequence[j].type)
            {
            case SEQUENCE_ITEM_VERTEX:
            case SEQUENCE_ITEM_FACE:
              BOOST_CHECK_EQUAL(endToStartCollector.m_sequence[j].index, startToEndCollector.m_sequence[k].index);
              break;
            case SEQUENCE_ITEM_EDGE:
              BOOST_CHECK_EQUAL(halfedgeIndexMap[endToStartCollector.m_sequence[j].halfedge], halfedgeIndexMap[CGAL::opposite(startToEndCollector.m_sequence[k].halfedge, P)]);
              break;
            }
          }
          else
          {
            Sequence_item<Traits>& steItem = startToEndCollector.m_sequence[k];
            Sequence_item<Traits>& etsItem = endToStartCollector.m_sequence[j];
            
            std::string names[2] = { "STE", "ETS" };
            Sequence_item<Traits> items[2] = { startToEndCollector.m_sequence[k], endToStartCollector.m_sequence[j] };
            
            for (size_t d = 0; d < 2; ++d)
            {
              if (items[d].type == SEQUENCE_ITEM_EDGE)
              {
                std::cout << "\t" << names[d] << "(edge): " << vertexIndexMap[CGAL::source(items[d].halfedge, P)] << " , " << vertexIndexMap[CGAL::target(items[d].halfedge, P)] << " : " << items[d].edgeAlpha << std::endl;
              }
              else if (items[d].type == SEQUENCE_ITEM_VERTEX)
              {
                std::cout << "\t" << names[d] << "(vertex): " << vertexIndexMap[items[d].vertex] << std::endl;
              }
            }
          }
        }
      }
    }

  }

  
}
  
// Hack to trick cgal_test_with_cmake into using this file even without a main
// int main(int argc, char** argv)
// {
//   return 0;
// }