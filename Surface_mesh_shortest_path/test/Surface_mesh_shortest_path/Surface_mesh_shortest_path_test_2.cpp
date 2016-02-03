#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <utility>
#include <cmath>

#include <CGAL/Random.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <CGAL/Surface_mesh_shortest_path/Surface_mesh_shortest_path_traits.h>
#include <CGAL/Surface_mesh_shortest_path/Surface_mesh_shortest_path.h>
#include <CGAL/Surface_mesh_shortest_path/function_objects.h>
#include <CGAL/Surface_mesh_shortest_path/barycentric.h>
#include <CGAL/Surface_mesh_shortest_path/internal/misc_functions.h>

#include <CGAL/test_util.h>
#include "check.h"

int main(int argc, char* argv[])
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;

  typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Polyhedron_3> Traits;
  typedef Traits::Barycentric_coordinate Barycentric_coordinate;
  typedef Traits::FT FT;
  typedef boost::graph_traits<Polyhedron_3> Graph_traits;
  typedef Graph_traits::vertex_descriptor vertex_descriptor;
  typedef Graph_traits::vertex_iterator vertex_iterator;
  typedef Graph_traits::face_descriptor face_descriptor;
  typedef Graph_traits::face_iterator face_iterator;
  typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
  typedef boost::property_map<Polyhedron_3, boost::vertex_index_t>::type VIM;
  typedef boost::property_map<Polyhedron_3, boost::halfedge_index_t>::type HIM;
  typedef boost::property_map<Polyhedron_3, boost::face_index_t>::type FIM;



  Traits traits;

  std::string mesh(argv[1]);

  int randSeed = 2681972;

  if (argc > 2)
  {
    randSeed = std::atoi(argv[2]);
  }

  CGAL::Random rand(randSeed);

  Polyhedron_3 polyhedron;
  std::ifstream in(mesh.c_str());

  in >> polyhedron;

  in.close();

  CGAL::set_halfedgeds_items_id(polyhedron);

  VIM vertexIndexMap(get(boost::vertex_index, polyhedron));
  HIM halfedgeIndexMap(get(boost::halfedge_index, polyhedron));
  FIM faceIndexMap(get(boost::face_index, polyhedron));

  vertex_iterator verticesStart;
  vertex_iterator verticesEnd;

  std::vector<vertex_descriptor> vertices;

  boost::tie(verticesStart, verticesEnd) = CGAL::vertices(polyhedron);

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
  Surface_mesh_shortest_path endToStartShortestPaths(polyhedron, traits);

  const size_t numTests = 15;

  std::cout << "Mesh: " << mesh << std::endl;


  for (size_t i = 0; i < numTests; ++i)
  {
    size_t startVertexIndex = rand.get_int(0, vertices.size());
    size_t endVertexIndex = rand.get_int(0, vertices.size());

    vertex_descriptor startVertex = vertices[startVertexIndex];
    vertex_descriptor endVertex = vertices[endVertexIndex];

    //startToEndShortestPaths.m_debugOutput = true;

    startToEndShortestPaths.clear();
    startToEndShortestPaths.add_source_point(startVertex);
    startToEndShortestPaths.build_sequence_tree();

    CGAL::test::Edge_sequence_collector<Traits> startToEndCollector(vertexIndexMap, halfedgeIndexMap, faceIndexMap);
    startToEndShortestPaths.shortest_path_sequence_to_source_points(endVertex, startToEndCollector);

    FT startToEnd = startToEndShortestPaths.shortest_distance_to_source_points(endVertex).first;

    //endToStartShortestPaths.m_debugOutput = true;

    endToStartShortestPaths.clear();
    endToStartShortestPaths.add_source_point(endVertex);
    endToStartShortestPaths.build_sequence_tree();

    CGAL::test::Edge_sequence_collector<Traits> endToStartCollector(vertexIndexMap, halfedgeIndexMap, faceIndexMap);
    endToStartShortestPaths.shortest_path_sequence_to_source_points(startVertex, endToStartCollector);

    FT endToStart = endToStartShortestPaths.shortest_distance_to_source_points(startVertex).first;

    CHECK_CLOSE(startToEnd, endToStart, FT(0.0000001));

    CHECK_EQUAL(startToEndCollector.m_sequence.size(), endToStartCollector.m_sequence.size());

    if (startToEndCollector.m_sequence.size() == endToStartCollector.m_sequence.size())
    {
      if (startToEndCollector.m_sequence.size() > 3)
      {
        size_t k = startToEndCollector.m_sequence.size() - 1;

        for (size_t j = 0; j < endToStartCollector.m_sequence.size(); ++j)
        {
          CHECK_EQUAL(endToStartCollector.m_sequence[j].type, startToEndCollector.m_sequence[k].type);

          if (endToStartCollector.m_sequence[j].type == startToEndCollector.m_sequence[k].type)
          {
            switch (endToStartCollector.m_sequence[j].type)
            {
            case CGAL::test::SEQUENCE_ITEM_VERTEX:
            case CGAL::test::SEQUENCE_ITEM_FACE:
              CHECK_EQUAL(endToStartCollector.m_sequence[j].index, startToEndCollector.m_sequence[k].index);
              break;
            case CGAL::test::SEQUENCE_ITEM_EDGE:
              CHECK_EQUAL(halfedgeIndexMap[endToStartCollector.m_sequence[j].halfedge], halfedgeIndexMap[CGAL::opposite(startToEndCollector.m_sequence[k].halfedge, polyhedron)]);
              break;
            }
          }
          else
          {
            std::string names[2] = { "STE", "ETS" };
            CGAL::test::Sequence_item<Traits> items[2] = { startToEndCollector.m_sequence[k], endToStartCollector.m_sequence[j] };
            Surface_mesh_shortest_path* pathStructures[2] = { &startToEndShortestPaths, &endToStartShortestPaths };

            for (size_t d = 0; d < 2; ++d)
            {
              if (items[d].type == CGAL::test::SEQUENCE_ITEM_EDGE)
              {
                std::cout << "\t" << names[d] << "(edge): " << vertexIndexMap[source(items[d].halfedge, polyhedron)] << " , " << vertexIndexMap[target(items[d].halfedge, polyhedron)] << " : " << items[d].edgeAlpha << std::endl;
              }
              else if (items[d].type == CGAL::test::SEQUENCE_ITEM_VERTEX)
              {
                std::cout << "\t" << names[d] << "(vertex): " << vertexIndexMap[items[d].vertex] << " , Distance: " << pathStructures[d]->shortest_distance_to_source_points(items[d].vertex).first << std::endl;
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

    Barycentric_coordinate startLocation = CGAL::test::random_coordinate<Traits>(rand);
    Barycentric_coordinate endLocation = CGAL::test::random_coordinate<Traits>(rand);

    //shortestPaths.m_debugOutput = true;

    startToEndShortestPaths.clear();
    startToEndShortestPaths.add_source_point(startFace, startLocation);
    startToEndShortestPaths.build_sequence_tree();

    FT startToEnd = startToEndShortestPaths.shortest_distance_to_source_points(endFace, endLocation).first;

    CGAL::test::Edge_sequence_collector<Traits> startToEndCollector(vertexIndexMap, halfedgeIndexMap, faceIndexMap);
    startToEndShortestPaths.shortest_path_sequence_to_source_points(endFace, endLocation, startToEndCollector);

    endToStartShortestPaths.clear();
    endToStartShortestPaths.add_source_point(endFace, endLocation);
    endToStartShortestPaths.build_sequence_tree();

    FT endToStart = endToStartShortestPaths.shortest_distance_to_source_points(startFace, startLocation).first;

    CGAL::test::Edge_sequence_collector<Traits> endToStartCollector(vertexIndexMap, halfedgeIndexMap, faceIndexMap);
    endToStartShortestPaths.shortest_path_sequence_to_source_points(startFace, startLocation, endToStartCollector);

    //std::cout << std::setprecision(15) << std::endl;

    CHECK_CLOSE(startToEnd, endToStart, FT(0.0000001));

    CHECK_EQUAL(startToEndCollector.m_sequence.size(), endToStartCollector.m_sequence.size());

    if (startToEndCollector.m_sequence.size() > 3 && startToEndCollector.m_sequence.size() == endToStartCollector.m_sequence.size())
    {

      for (size_t j = 0; j < endToStartCollector.m_sequence.size(); ++j)
      {
        size_t k = endToStartCollector.m_sequence.size() - j - 1;

        CHECK_EQUAL(endToStartCollector.m_sequence[j].type, startToEndCollector.m_sequence[k].type);

        if (endToStartCollector.m_sequence[j].type == startToEndCollector.m_sequence[k].type)
        {
          switch (endToStartCollector.m_sequence[j].type)
          {
          case CGAL::test::SEQUENCE_ITEM_VERTEX:
          case CGAL::test::SEQUENCE_ITEM_FACE:
            CHECK_EQUAL(endToStartCollector.m_sequence[j].index, startToEndCollector.m_sequence[k].index);
            break;
          case CGAL::test::SEQUENCE_ITEM_EDGE:
            CHECK_EQUAL(halfedgeIndexMap[endToStartCollector.m_sequence[j].halfedge], halfedgeIndexMap[CGAL::opposite(startToEndCollector.m_sequence[k].halfedge, polyhedron)]);
            break;
          }
        }
        else
        {
          std::string names[2] = { "STE", "ETS" };
          CGAL::test::Sequence_item<Traits> items[2] = { startToEndCollector.m_sequence[k], endToStartCollector.m_sequence[j] };

          for (size_t d = 0; d < 2; ++d)
          {
            if (items[d].type == CGAL::test::SEQUENCE_ITEM_EDGE)
            {
              std::cout << "\t" << names[d] << "(edge): " << vertexIndexMap[CGAL::source(items[d].halfedge, polyhedron)] << " , " << vertexIndexMap[CGAL::target(items[d].halfedge, polyhedron)] << " : " << items[d].edgeAlpha << std::endl;
            }
            else if (items[d].type == CGAL::test::SEQUENCE_ITEM_VERTEX)
            {
              std::cout << "\t" << names[d] << "(vertex): " << vertexIndexMap[items[d].vertex] << std::endl;
            }
          }
        }
      }
    }
  }

  return 0;
}


