#include <iostream>
#include <fstream>
#include <utility>
#include <set>

#include <CGAL/Random.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/boost/graph/iterator.h>

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
  typedef Traits::FT FT;
  typedef boost::graph_traits<Polyhedron_3> Graph_traits;
  typedef Graph_traits::face_descriptor face_descriptor;
  typedef Graph_traits::face_iterator face_iterator;
  typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
  typedef Surface_mesh_shortest_path::Face_location Face_location;
  typedef boost::property_map<Polyhedron_3, boost::vertex_index_t>::const_type VIM;
  typedef boost::property_map<Polyhedron_3, boost::halfedge_index_t>::const_type HIM;
  typedef boost::property_map<Polyhedron_3, boost::face_index_t>::const_type FIM;

  Traits traits;

  std::string mesh(argv[1]);

  int randSeed = 4983304;

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

  CGAL_USE(vertexIndexMap);
  CGAL_USE(halfedgeIndexMap);
  CGAL_USE(faceIndexMap);

  face_iterator facesStart;
  face_iterator facesEnd;

  std::vector<face_descriptor> faces;

  boost::tie(facesStart, facesEnd) = CGAL::faces(polyhedron);

  for (face_iterator it = facesStart; it != facesEnd; ++it)
  {
    faces.push_back(*it);
  }

  Surface_mesh_shortest_path shortestPaths(polyhedron, traits);

  const size_t numInitialLocations = 10;

  std::vector<Face_location> sourcePoints;
  std::vector<Surface_mesh_shortest_path::Source_point_iterator> handles;

  // First, try adding a few locations

  for (size_t i = 0; i < numInitialLocations; ++i)
  {
    size_t faceId = static_cast<std::size_t>(rand.get_int(0, static_cast<int>(faces.size())));
    sourcePoints.push_back(Face_location(faces[faceId], CGAL::test::random_coordinates<Traits>(rand)));
    shortestPaths.add_source_point(sourcePoints.back().first, sourcePoints.back().second);
  }

  CHECK_EQUAL(numInitialLocations, shortestPaths.number_of_source_points());

  size_t checkNumLocations = 0;

  for (Surface_mesh_shortest_path::Source_point_iterator it = shortestPaths.source_points_begin(); it != shortestPaths.source_points_end(); ++it)
  {
    handles.push_back(it);
    ++checkNumLocations;
  }

  CHECK_EQUAL(checkNumLocations, shortestPaths.number_of_source_points());

  for (Surface_mesh_shortest_path::Source_point_iterator it = shortestPaths.source_points_begin(); it != shortestPaths.source_points_end(); ++it)
  {
    Surface_mesh_shortest_path::Shortest_path_result result = shortestPaths.shortest_distance_to_source_points(it->first, it->second);

    CHECK_CLOSE(FT(0.0), result.first, FT(0.000001));
    assert(result.second == it);
  }

  // Then, remove half of them

  for (size_t i = 0; i < handles.size(); ++i)
  {
    if (i % 2 == 0)
    {
      shortestPaths.remove_source_point(handles[i]);
    }
  }

  CHECK_EQUAL(numInitialLocations / 2, shortestPaths.number_of_source_points());

  // and ensure that they are indeed removed
  for (size_t i = 0; i < sourcePoints.size(); ++i)
  {
    Surface_mesh_shortest_path::Shortest_path_result result = shortestPaths.shortest_distance_to_source_points(sourcePoints[i].first, sourcePoints[i].second);

    if (i % 2 != 0)
    {
      CHECK_CLOSE(FT(0.0), result.first, FT(0.000001));
      assert(handles[i] == result.second);
    }
    else
    {
      if ( !(result.first < FT(0.0) || result.first > FT(0.00001) ) )
      {
        std::cerr << "Incorrect resulting distance: " << result.first << "\n";
        return 1;
      }
    }
  }

  // add a few back
  for (size_t i = 0; i < handles.size(); ++i)
  {
    if (i % 2 == 0)
    {
      handles[i] = shortestPaths.add_source_point(sourcePoints[i]);
    }
  }

  // ... and remove some others
  for (size_t i = 0; i < handles.size(); ++i)
  {
    if (i % 3 == 0)
    {
      shortestPaths.remove_source_point(handles[i]);
    }
  }

  // and check it once again
  for (size_t i = 0; i < sourcePoints.size(); ++i)
  {
    Surface_mesh_shortest_path::Shortest_path_result result = shortestPaths.shortest_distance_to_source_points(sourcePoints[i].first, sourcePoints[i].second);

    if (i % 3 != 0)
    {
      CHECK_CLOSE(FT(0.0), result.first, FT(0.000001));
      assert(handles[i] == result.second);
    }
    else
    {
      if( !(result.first < FT(0.0) || result.first > FT(0.00001)) )
      {
        std::cerr << "Resulted distance: " << result.first << "\n";
        return 1;
      }
    }
  }
   return 0;
}



