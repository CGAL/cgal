#include <iostream>
#include <fstream>
#include <utility>

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

#include <CGAL/use.h>

#include <CGAL/test_util.h>
#include "check.h"

void shortest_path_regular_tetrahedron()
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;
  typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Polyhedron_3> Traits;
  typedef Traits::Barycentric_coordinates Barycentric_coordinates;
  typedef Traits::FT FT;
  typedef boost::graph_traits<Polyhedron_3> Graph_traits;
  typedef Graph_traits::vertex_iterator vertex_iterator;
  typedef Graph_traits::face_descriptor face_descriptor;
  typedef Graph_traits::face_iterator face_iterator;
  typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;

  Traits traits;

  Traits::Construct_barycentric_coordinates construct_barycentric_coordinates(traits.construct_barycentric_coordinates_object());

  Polyhedron_3 P;

  CGAL::test::make_regular_tetrahedron(P);

  CGAL::set_halfedgeds_items_id(P);

  Barycentric_coordinates b = construct_barycentric_coordinates(FT(1.0) / FT(3.0), FT(1.0) / FT(3.0), FT(1.0) / FT(3.0));

  face_iterator startFace;
  face_iterator endFace;

  boost::tie(startFace,endFace) = CGAL::faces(P);

  face_descriptor firstFace = *startFace;

  Surface_mesh_shortest_path shortestPaths(P, traits);
  //shortestPaths.m_debugOutput = true;
  shortestPaths.add_source_point(firstFace, b);
  shortestPaths.build_sequence_tree();

  vertex_iterator currentVertex;
  vertex_iterator endVertex;

  Kernel::FT sideLength = Kernel::FT(2.0);
  Kernel::FT halfSideLength = sideLength / Kernel::FT(2.0);
  Kernel::FT triangleHeight = CGAL::sqrt((sideLength*sideLength) - (halfSideLength*halfSideLength));

  for (boost::tie(currentVertex, endVertex) = CGAL::vertices(P); currentVertex != endVertex; ++currentVertex)
  {
    if ((*currentVertex)->point().y()==-1)
    {
      CHECK_CLOSE(shortestPaths.shortest_distance_to_source_points(*currentVertex).first, Kernel::FT((triangleHeight * Kernel::FT(4.0)) / Kernel::FT(3.0)), Kernel::FT(0.000001));
    }
    else
    {
      CHECK_CLOSE(shortestPaths.shortest_distance_to_source_points(*currentVertex).first, Kernel::FT((triangleHeight * Kernel::FT(2.0)) / Kernel::FT(3.0)), Kernel::FT(0.000001));
    }
  }
}

void test_simple_saddle_vertex_mesh()
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;
  typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Polyhedron_3> Traits;
  typedef Traits::Barycentric_coordinates Barycentric_coordinates;
  typedef Traits::FT FT;
  typedef Traits::Point_3 Point_3;
  typedef Traits::Point_2 Point_2;
  typedef Traits::Triangle_3 Triangle_3;
  typedef Traits::Triangle_2 Triangle_2;
  typedef Traits::Segment_2 Segment_2;
  typedef boost::graph_traits<Polyhedron_3> Graph_traits;
  typedef Graph_traits::vertex_descriptor vertex_descriptor;
  typedef Graph_traits::vertex_iterator vertex_iterator;
  typedef Graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef Graph_traits::face_descriptor face_descriptor;
  typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
  typedef boost::property_map<Polyhedron_3, CGAL::vertex_point_t>::const_type VPM;

  Traits traits;

  Traits::Compute_squared_distance_3 compute_squared_distance_3(traits.compute_squared_distance_3_object());
  Traits::Compute_squared_distance_2 compute_squared_distance_2(traits.compute_squared_distance_2_object());
  Traits::Construct_triangle_3_along_segment_2_flattening flatten_triangle_3_along_segment_2(traits.construct_triangle_3_along_segment_2_flattening_object());
  Traits::Construct_barycentric_coordinates construct_barycentric_coordinates(traits.construct_barycentric_coordinates_object());

  std::ifstream inFile("data/saddle_vertex_mesh.off");

  Polyhedron_3 P;

  inFile >> P;

  inFile.close();

  CGAL::set_halfedgeds_items_id(P);

  vertex_iterator startVertex;
  vertex_iterator endVertex;
  boost::tie(startVertex, endVertex) = CGAL::vertices(P);

  vertex_iterator currentVertex = startVertex;

  ++currentVertex;
  vertex_descriptor rootSearchVertex = *currentVertex;

  face_descriptor currentFace = CGAL::face(CGAL::halfedge(rootSearchVertex, P), P);
  size_t vertexIndex = CGAL::test::face_vertex_index(currentFace, rootSearchVertex, P);
  Barycentric_coordinates baryCoord = construct_barycentric_coordinates(vertexIndex == 0 ? FT(1.0) : FT(0.0), vertexIndex == 1 ? FT(1.0) : FT(0.0), vertexIndex == 2 ? FT(1.0) : FT(0.0));

  Surface_mesh_shortest_path shortestPaths(P, traits);
  //shortestPaths.m_debugOutput = true;
  Surface_mesh_shortest_path::Source_point_iterator firstSourcePoint = shortestPaths.add_source_point(currentFace, baryCoord);
  shortestPaths.build_sequence_tree();

  VPM vpm = CGAL::get(CGAL::vertex_point, P);

  Point_3 vertexLocations[8];
  vertex_descriptor vertexHandles[8];

  currentVertex = startVertex;

  for (size_t i = 0; i < 8; ++i)
  {
    vertexHandles[i] = *currentVertex;
    vertexLocations[i] = vpm[*currentVertex];
    ++currentVertex;
  }

  FT distanceToBottom = CGAL::sqrt(compute_squared_distance_3(vertexLocations[1], vertexLocations[0]));
  FT largerSideLength = CGAL::sqrt(compute_squared_distance_3(vertexLocations[1], vertexLocations[2]));
  FT distanceToSaddle = CGAL::sqrt(compute_squared_distance_3(vertexLocations[1], vertexLocations[4]));
  FT shorterSideLength = CGAL::sqrt(compute_squared_distance_3(vertexLocations[4], vertexLocations[6]));

  Triangle_3 lower(vertexLocations[6], vertexLocations[4], vertexLocations[1]);
  Triangle_3 upper(vertexLocations[4], vertexLocations[6], vertexLocations[7]);

  Segment_2 base(Point_2(CGAL::ORIGIN), Point_2(shorterSideLength, FT(0.0)));

  Triangle_2 flatLower(flatten_triangle_3_along_segment_2(lower, 0, base));
  Triangle_2 flatUpper(flatten_triangle_3_along_segment_2(upper, 0, Segment_2(base[1], base[0])));

  FT distanceToApex = CGAL::sqrt(compute_squared_distance_2(flatLower[2], flatUpper[2]));

  FT expectedDistances[8] =
  {
    distanceToBottom, // a vertex of the larger tetrahedron
    FT(0.0), // the initial vertex
    largerSideLength, // a vertex of the larger tetrahedron
    largerSideLength, // a vertex of the larger tetrahedron
    distanceToSaddle,  // direct line of sight from root
    distanceToSaddle + shorterSideLength,  // around the corner from a pseudo-source (not in direct line of geodesic sight)
    distanceToSaddle, // direct line of sight from root
    distanceToApex,
  };

  currentVertex = startVertex;

  for (size_t i = 0; i < 8; ++i)
  {
    CHECK_CLOSE(shortestPaths.shortest_distance_to_source_points(*currentVertex).first, expectedDistances[i], Kernel::FT(0.0001));
    ++currentVertex;
  }

  // test the edge sequence reporting
  CGAL::test::Edge_sequence_collector<Traits> collector(P);

  shortestPaths.shortest_path_sequence_to_source_points(vertexHandles[5], collector);

  CHECK_EQUAL(collector.m_sequence.size(), 3u);
  CHECK_EQUAL(collector.m_sequence[1].type, CGAL::test::SEQUENCE_ITEM_VERTEX);
  assert(collector.m_sequence[1].index == 4 || collector.m_sequence[1].index == 6);
  CHECK_EQUAL(collector.m_sequence[2].type, CGAL::test::SEQUENCE_ITEM_VERTEX);
  CHECK_EQUAL(collector.m_sequence[2].index, 1u);

  collector.m_sequence.clear();

  typedef boost::property_map<Polyhedron_3, CGAL::halfedge_external_index_t>::const_type HalfedgeIndexMap;

  HalfedgeIndexMap halfedgeIndexMap(CGAL::get(CGAL::halfedge_external_index, P));

  shortestPaths.shortest_path_sequence_to_source_points(vertexHandles[7], collector);

  CHECK_EQUAL(collector.m_sequence.size(), 3u);
  CHECK_EQUAL(collector.m_sequence[1].type, CGAL::test::SEQUENCE_ITEM_EDGE);
  CHECK_EQUAL(collector.m_sequence[1].index, halfedgeIndexMap[CGAL::halfedge(vertexHandles[4], vertexHandles[6], P).first]);
  CHECK_CLOSE(collector.m_sequence[1].edgeAlpha, FT(0.5), FT(0.0001));
  CHECK_EQUAL(collector.m_sequence[2].type, CGAL::test::SEQUENCE_ITEM_VERTEX);
  CHECK_EQUAL(collector.m_sequence[2].index, 1u);

  // Now test an internal face location sequence
  halfedge_descriptor firstCrossing = CGAL::halfedge(vertexHandles[4], vertexHandles[7], P).first;

  size_t edgeIndex = CGAL::Surface_mesh_shortest_paths_3::internal::edge_index(firstCrossing, P);

  Barycentric_coordinates location = construct_barycentric_coordinates(0.25, 0.5, 0.25);

  collector.m_sequence.clear();
  shortestPaths.shortest_path_sequence_to_source_points(CGAL::face(firstCrossing, P), construct_barycentric_coordinates(location[edgeIndex], location[(edgeIndex + 1) % 3], location[(edgeIndex + 2) % 3]), collector);

  CHECK_EQUAL(collector.m_sequence.size(), 4u);
  CHECK_EQUAL(collector.m_sequence[1].type, CGAL::test::SEQUENCE_ITEM_EDGE);
  CHECK_EQUAL(collector.m_sequence[1].index, halfedgeIndexMap[firstCrossing]);
  CHECK_EQUAL(collector.m_sequence[2].type, CGAL::test::SEQUENCE_ITEM_EDGE);
  CHECK_EQUAL(collector.m_sequence[2].index, halfedgeIndexMap[CGAL::halfedge(vertexHandles[4], vertexHandles[6], P).first]);
  CHECK_EQUAL(collector.m_sequence[3].type, CGAL::test::SEQUENCE_ITEM_VERTEX);
  CHECK_EQUAL(collector.m_sequence[3].index, 1u);

  // Now test with 2 source vertices
  currentVertex = startVertex;

  for (size_t i = 0; i < 5; ++i)
  {
    ++currentVertex;
  }

  vertex_descriptor rootSearchVertex2 = *currentVertex;

  face_descriptor currentFace2 = CGAL::face(CGAL::halfedge(rootSearchVertex2, P), P);
  size_t vertexIndex2 = CGAL::test::face_vertex_index(currentFace2, rootSearchVertex2, P);
  Barycentric_coordinates baryCoord2 = construct_barycentric_coordinates(vertexIndex2 == 0 ? FT(1.0) : FT(0.0), vertexIndex2 == 1 ? FT(1.0) : FT(0.0), vertexIndex2 == 2 ? FT(1.0) : FT(0.0));

  Surface_mesh_shortest_path::Source_point_iterator secondSourcePoint = shortestPaths.add_source_point(currentFace2, baryCoord2);
  shortestPaths.build_sequence_tree();

  FT distanceToApexFrom2 = CGAL::sqrt(compute_squared_distance_3(vertexLocations[5], vertexLocations[7]));

  Triangle_3 lower2(vertexLocations[2], vertexLocations[3], vertexLocations[5]);
  Triangle_3 upper2(vertexLocations[3], vertexLocations[2], vertexLocations[0]);

  Segment_2 base2(Point_2(CGAL::ORIGIN), Point_2(largerSideLength, FT(0.0)));

  Triangle_2 flatLower2(flatten_triangle_3_along_segment_2(lower2, 0, base2));
  Triangle_2 flatUpper2(flatten_triangle_3_along_segment_2(upper2, 0, Segment_2(base2[1], base2[0])));

  FT distanceToBottom2 = CGAL::sqrt(compute_squared_distance_2(flatLower2[2], flatUpper2[2]));

  FT expectedDistances2[8] =
  {
    distanceToBottom2, // a vertex of the larger tetrahedron
    FT(0.0), // an initial vertex
    distanceToSaddle, // a vertex of the larger tetrahedron
    distanceToSaddle, // a vertex of the larger tetrahedron
    shorterSideLength,  // direct line of sight from root
    FT(0.0),  // around the corner from a pseudo-source (not in direct line of geodesic sight)
    shorterSideLength, // direct line of sight from root
    distanceToApexFrom2,
  };

  Surface_mesh_shortest_path::Source_point_iterator expectedSources2[8] =
  {
    secondSourcePoint,
    firstSourcePoint,
    secondSourcePoint,
    secondSourcePoint,
    secondSourcePoint,
    secondSourcePoint,
    secondSourcePoint,
    secondSourcePoint,
  };

  currentVertex = startVertex;

  for (size_t i = 0; i < 8; ++i)
  {
    Surface_mesh_shortest_path::Shortest_path_result result = shortestPaths.shortest_distance_to_source_points(*currentVertex);

    CHECK_CLOSE(result.first, expectedDistances2[i], Kernel::FT(0.0001));
    assert(result.second == expectedSources2[i]);
    CGAL_USE(expectedSources2);
    ++currentVertex;
  }

  // Test removing a source vertex

  shortestPaths.remove_source_point(firstSourcePoint);
  shortestPaths.build_sequence_tree();

  // replace the only shortest path which was not reporting to the 2nd source point
  expectedDistances2[1] = distanceToSaddle + shorterSideLength;

  currentVertex = startVertex;

  for (size_t i = 0; i < 8; ++i)
  {
    Surface_mesh_shortest_path::Shortest_path_result result = shortestPaths.shortest_distance_to_source_points(*currentVertex);

    CHECK_CLOSE(result.first, expectedDistances2[i], Kernel::FT(0.0001));
    assert(result.second == secondSourcePoint);
    ++currentVertex;
  }

}

void test_boundary_mesh()
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;
  typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Polyhedron_3> Traits;
  typedef Traits::Barycentric_coordinates Barycentric_coordinates;
  typedef Traits::FT FT;
  typedef Traits::Point_3 Point_3;
  typedef Traits::Triangle_3 Triangle_3;
  typedef boost::graph_traits<Polyhedron_3> Graph_traits;
  typedef Graph_traits::vertex_descriptor vertex_descriptor;
  typedef Graph_traits::vertex_iterator vertex_iterator;
  typedef Graph_traits::face_descriptor face_descriptor;
  typedef Graph_traits::face_iterator face_iterator;
  typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
  typedef boost::property_map<Polyhedron_3, CGAL::vertex_point_t>::const_type VPM;

  Traits traits;

  Traits::Construct_triangle_3_to_triangle_2_projection project_triangle_3_to_triangle_2(traits.construct_triangle_3_to_triangle_2_projection_object());
  CGAL_USE(project_triangle_3_to_triangle_2);
  Traits::Compute_squared_distance_3 compute_squared_distance_3(traits.compute_squared_distance_3_object());
  Traits::Construct_barycenter_3 construct_barycenter_3(traits.construct_barycenter_3_object());
  Traits::Construct_triangle_3_along_segment_2_flattening flatten_triangle_3_along_segment_2(traits.construct_triangle_3_along_segment_2_flattening_object());
  CGAL_USE(flatten_triangle_3_along_segment_2);
  Traits::Construct_barycentric_coordinates construct_barycentric_coordinates(traits.construct_barycentric_coordinates_object());

  struct Construct_barycenter_in_triangle_3
  {
    Traits::Construct_barycenter_3 m_cb3;

    Construct_barycenter_in_triangle_3(Traits::Construct_barycenter_3 cb3)
      : m_cb3(cb3)
    {
    }

    Point_3 operator() (const Triangle_3& t, const Barycentric_coordinates& b)
    {
      return m_cb3(t[0], b[0], t[1], b[1], t[2], b[2]);
    }
  } construct_barycenter_in_triangle_3(construct_barycenter_3);

  std::ifstream inFile("data/boundary_mesh.off");

  Polyhedron_3 P;

  inFile >> P;

  inFile.close();

  CGAL::set_halfedgeds_items_id(P);

  face_iterator startFace;
  face_iterator endFace;

  boost::tie(startFace, endFace) = CGAL::faces(P);

  vertex_iterator currentVertex;
  vertex_iterator endVertex;

  VPM vpm = CGAL::get(CGAL::vertex_point, P);

  vertex_descriptor vertexHandles[10];
  face_descriptor faceHandles[8];
  Point_3 vertexLocations[10];
  size_t currentVertexIndex = 0;

  for (boost::tie(currentVertex, endVertex) = CGAL::vertices(P); currentVertex != endVertex; ++currentVertex)
  {
    vertexHandles[currentVertexIndex] = *currentVertex;
    vertexLocations[currentVertexIndex] = vpm[*currentVertex];
    ++currentVertexIndex;
  }

  size_t currentFaceIndex = 0;

  for (face_iterator currentFace = startFace; currentFace != endFace; ++currentFace)
  {
    faceHandles[currentFaceIndex] = *currentFace;
    ++currentFaceIndex;
  }

  Barycentric_coordinates startLocation = construct_barycentric_coordinates(FT(0.1), FT(0.8), FT(0.1));

  typedef boost::property_map<Polyhedron_3, CGAL::face_external_index_t>::const_type FaceIndexMap;

  FaceIndexMap faceIndexMap(CGAL::get(CGAL::face_external_index, P));

  Surface_mesh_shortest_path shortestPaths(P, traits);
  //shortestPaths.m_debugOutput = true;
  shortestPaths.add_source_point(*startFace, startLocation);
  shortestPaths.build_sequence_tree();

  Triangle_3 firstTriangle(vertexLocations[1], vertexLocations[0], vertexLocations[2]);

  Point_3 locationInTriangle(construct_barycenter_in_triangle_3(firstTriangle, startLocation));

  FT dist0 = shortestPaths.shortest_distance_to_source_points(vertexHandles[0]).first;
  CHECK_CLOSE(dist0, CGAL::sqrt(compute_squared_distance_3(locationInTriangle, vertexLocations[0])), FT(0.000001));

  FT dist1 = shortestPaths.shortest_distance_to_source_points(vertexHandles[1]).first;
  CHECK_CLOSE(dist1, CGAL::sqrt(compute_squared_distance_3(locationInTriangle, vertexLocations[1])), FT(0.000001));

  FT dist2 = shortestPaths.shortest_distance_to_source_points(vertexHandles[2]).first;
  CHECK_CLOSE(dist2, CGAL::sqrt(compute_squared_distance_3(locationInTriangle, vertexLocations[2])), FT(0.000001));

  FT dist3 = shortestPaths.shortest_distance_to_source_points(vertexHandles[3]).first;
  CHECK_CLOSE(dist3, CGAL::sqrt(compute_squared_distance_3(locationInTriangle, vertexLocations[3])), FT(0.000001));

  FT dist4 = shortestPaths.shortest_distance_to_source_points(vertexHandles[4]).first;
  CHECK_CLOSE(dist4, CGAL::sqrt(compute_squared_distance_3(locationInTriangle, vertexLocations[1])) + CGAL::sqrt(compute_squared_distance_3(vertexLocations[1], vertexLocations[4])), FT(0.000001));

  FT dist5 = shortestPaths.shortest_distance_to_source_points(vertexHandles[5]).first;
  CHECK_CLOSE(dist5, CGAL::sqrt(compute_squared_distance_3(locationInTriangle, vertexLocations[3])) + CGAL::sqrt(compute_squared_distance_3(vertexLocations[3], vertexLocations[5])), FT(0.000001));

  Barycentric_coordinates somewhereElseInFirstTriangle = construct_barycentric_coordinates(0.8, 0.05, 0.15);

  FT distT0 = shortestPaths.shortest_distance_to_source_points(faceHandles[0], somewhereElseInFirstTriangle).first;
  CHECK_CLOSE(distT0, CGAL::sqrt(compute_squared_distance_3(locationInTriangle, construct_barycenter_in_triangle_3(firstTriangle, somewhereElseInFirstTriangle))), FT(0.000001));

  Triangle_3 oneStepTriangle(vertexLocations[4], vertexLocations[1], vertexLocations[3]);
  Barycentric_coordinates locationInOneStepTriangle = construct_barycentric_coordinates(0.1, 0.8, 0.1);

  CGAL::test::Edge_sequence_collector<Traits> collector(P);
  shortestPaths.shortest_path_sequence_to_source_points(faceHandles[2], locationInOneStepTriangle, collector);

  FT distT2 = shortestPaths.shortest_distance_to_source_points(faceHandles[2], locationInOneStepTriangle).first;
  CHECK_CLOSE(distT2, dist1 + CGAL::sqrt(compute_squared_distance_3(vertexLocations[1], construct_barycenter_in_triangle_3(oneStepTriangle, locationInOneStepTriangle))), FT(0.00001));

  Triangle_3 twoStepTriangle(vertexLocations[6], vertexLocations[5], vertexLocations[7]);
  Barycentric_coordinates locationInTwoStepTriangle = construct_barycentric_coordinates(0.8, 0.1, 0.1);

  FT distT5 = shortestPaths.shortest_distance_to_source_points(faceHandles[5], locationInTwoStepTriangle).first;
  CHECK_CLOSE(distT5, dist3 + CGAL::sqrt(compute_squared_distance_3(vertexLocations[3], construct_barycenter_in_triangle_3(twoStepTriangle, locationInTwoStepTriangle))), FT(0.00001));

  Triangle_3 threeStepTriangle(vertexLocations[7], vertexLocations[5], vertexLocations[8]);
  Barycentric_coordinates locationInThreeStepTriangle = construct_barycentric_coordinates(0.2, 0.6, 0.2);

  FT distT6 = shortestPaths.shortest_distance_to_source_points(faceHandles[6], locationInThreeStepTriangle).first;
  CHECK_CLOSE(distT6, dist5 + CGAL::sqrt(compute_squared_distance_3(vertexLocations[5], construct_barycenter_in_triangle_3(threeStepTriangle, locationInThreeStepTriangle))), FT(0.00001));
}


int main()
{
  shortest_path_regular_tetrahedron();
  test_simple_saddle_vertex_mesh();
  test_boundary_mesh();

   return 0;
}
