// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#include <CGAL/Interval_nt.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

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

#define BOOST_TEST_MODULE polyhedron_shortest_path_test
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE( shortest_path_regular_tetrahedron )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
  typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Polyhedron_3> Traits;
  typedef Traits::Barycentric_coordinate Barycentric_coordinate;
  typedef Traits::FT FT;
  typedef boost::graph_traits<Polyhedron_3> GraphTraits;
  typedef GraphTraits::vertex_descriptor vertex_descriptor;
  typedef GraphTraits::vertex_iterator vertex_iterator;
  typedef GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef GraphTraits::halfedge_iterator halfedge_iterator;
  typedef GraphTraits::face_descriptor face_descriptor;
  typedef GraphTraits::face_iterator face_iterator;
  typedef CGAL::Polyhedron_shortest_path<Traits> Polyhedron_shortest_path;
  
  Traits traits;

  Polyhedron_3 P;
  
  CGAL::test::make_regular_tetrahedron(P);
  
  Barycentric_coordinate b(FT(1.0) / FT(3.0), FT(1.0) / FT(3.0), FT(1.0) / FT(3.0));
  
  face_iterator startFace;
  face_iterator endFace;
  
  boost::tie(startFace,endFace) = CGAL::faces(P);
  
  face_descriptor firstFace = *startFace;
  
  Polyhedron_shortest_path shortestPaths(P, traits);
  //shortestPaths.m_debugOutput = true;
  shortestPaths.compute_shortest_paths(firstFace, b);

  vertex_iterator currentVertex;
  vertex_iterator endVertex;
  
  size_t vertexIndex = 0;
  
  Kernel::FT sideLength = Kernel::FT(2.0);
  Kernel::FT halfSideLength = sideLength / Kernel::FT(2.0);
  Kernel::FT triangleHeight = CGAL::sqrt((sideLength*sideLength) - (halfSideLength*halfSideLength));

  for (boost::tie(currentVertex, endVertex) = CGAL::vertices(P); currentVertex != endVertex; ++currentVertex)
  {
    if (vertexIndex == 0)
    {
      BOOST_CHECK_CLOSE(shortestPaths.shortest_distance_to_vertex(*currentVertex).first, Kernel::FT((triangleHeight * Kernel::FT(4.0)) / Kernel::FT(3.0)), Kernel::FT(0.000001));
    }
    else
    {
      BOOST_CHECK_CLOSE(shortestPaths.shortest_distance_to_vertex(*currentVertex).first, Kernel::FT((triangleHeight * Kernel::FT(2.0)) / Kernel::FT(3.0)), Kernel::FT(0.000001));
    }
    
    ++vertexIndex;
  }
}
  
BOOST_AUTO_TEST_CASE( test_simple_saddle_vertex_mesh )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
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

  Traits traits;

  Traits::Compute_squared_distance_3 compute_squared_distance_3(traits.compute_squared_distance_3_object());
  Traits::Compute_squared_distance_2 compute_squared_distance_2(traits.compute_squared_distance_2_object());
  Traits::Flatten_triangle_3_along_segment_2 flatten_triangle_3_along_segment_2(traits.flatten_triangle_3_along_segment_2_object());
  
  std::ifstream inFile("data/saddle_vertex_mesh.off");
  
  Polyhedron_3 P;
  
  inFile >> P;

  inFile.close();
  
  vertex_iterator startVertex;
  vertex_iterator endVertex;
  boost::tie(startVertex, endVertex) = CGAL::vertices(P);
  
  vertex_iterator currentVertex = startVertex;
  
  ++currentVertex;
  vertex_descriptor rootSearchVertex = *currentVertex;
  
  face_descriptor currentFace = CGAL::face(CGAL::halfedge(rootSearchVertex, P), P);
  size_t vertexIndex = CGAL::test::face_vertex_index(currentFace, rootSearchVertex, P);
  Barycentric_coordinate baryCoord(vertexIndex == 0 ? FT(1.0) : FT(0.0), vertexIndex == 1 ? FT(1.0) : FT(0.0), vertexIndex == 2 ? FT(1.0) : FT(0.0));

  Polyhedron_shortest_path shortestPaths(P, traits);
  //shortestPaths.m_debugOutput = true;
  shortestPaths.compute_shortest_paths(currentFace, baryCoord);

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
    BOOST_CHECK_CLOSE(shortestPaths.shortest_distance_to_vertex(*currentVertex).first, expectedDistances[i], Kernel::FT(0.0001));
    ++currentVertex;
  }
  
  // test the edge sequence reporting
  CGAL::test::Edge_sequence_collector<Traits> collector(P);
  
  shortestPaths.shortest_path_sequence(vertexHandles[5], collector);
  
  BOOST_CHECK_EQUAL(collector.m_sequence.size(), 2);
  BOOST_CHECK_EQUAL(collector.m_sequence[0].type, CGAL::test::SEQUENCE_ITEM_VERTEX);
  //BOOST_CHECK_EQUAL(collector.m_sequence[0].index, 4 || collector.m_sequence[0].index, 6);
  BOOST_CHECK_EQUAL(collector.m_sequence[1].type, CGAL::test::SEQUENCE_ITEM_VERTEX);
  BOOST_CHECK_EQUAL(collector.m_sequence[1].index, 1);
  
  collector.m_sequence.clear();
  
  typedef boost::property_map<Polyhedron_3, CGAL::halfedge_external_index_t>::type HalfedgeIndexMap;
  
  HalfedgeIndexMap halfedgeIndexMap(CGAL::get(CGAL::halfedge_external_index, P));
  
  shortestPaths.shortest_path_sequence(vertexHandles[7], collector);
  
  BOOST_CHECK_EQUAL(collector.m_sequence.size(), 2);
  BOOST_CHECK_EQUAL(collector.m_sequence[0].type, CGAL::test::SEQUENCE_ITEM_EDGE);
  BOOST_CHECK_EQUAL(collector.m_sequence[0].index, halfedgeIndexMap[CGAL::halfedge(vertexHandles[4], vertexHandles[6], P).first]);
  BOOST_CHECK_CLOSE(collector.m_sequence[0].edgeAlpha, FT(0.5), FT(0.0001));
  BOOST_CHECK_EQUAL(collector.m_sequence[1].type, CGAL::test::SEQUENCE_ITEM_VERTEX);
  BOOST_CHECK_EQUAL(collector.m_sequence[1].index, 1);
  
  // Now test an internal face location sequence
  halfedge_descriptor firstCrossing = CGAL::halfedge(vertexHandles[4], vertexHandles[7], P).first;
  
  size_t edgeIndex = CGAL::internal::edge_index(firstCrossing, P);

  Barycentric_coordinate location(0.25, 0.5, 0.25);

  collector.m_sequence.clear();
  shortestPaths.shortest_path_sequence(CGAL::face(firstCrossing, P), Traits::Barycentric_coordinate(location[edgeIndex], location[(edgeIndex + 1) % 3], location[(edgeIndex + 2) % 3]), collector);

  BOOST_CHECK_EQUAL(collector.m_sequence.size(), 3);
  BOOST_CHECK_EQUAL(collector.m_sequence[0].type, CGAL::test::SEQUENCE_ITEM_EDGE);
  BOOST_CHECK_EQUAL(collector.m_sequence[0].index, halfedgeIndexMap[firstCrossing]);
  BOOST_CHECK_EQUAL(collector.m_sequence[1].type, CGAL::test::SEQUENCE_ITEM_EDGE);
  BOOST_CHECK_EQUAL(collector.m_sequence[1].index, halfedgeIndexMap[CGAL::halfedge(vertexHandles[4], vertexHandles[6], P).first]);
  BOOST_CHECK_EQUAL(collector.m_sequence[2].type, CGAL::test::SEQUENCE_ITEM_VERTEX);
  BOOST_CHECK_EQUAL(collector.m_sequence[2].index, 1);
  
  // Now test with 2 source vertices
  currentVertex = startVertex;
  
  for (size_t i = 0; i < 5; ++i)
  {
    ++currentVertex;
  }
  
  vertex_descriptor rootSearchVertex2 = *currentVertex;
  
  face_descriptor currentFace2 = CGAL::face(CGAL::halfedge(rootSearchVertex2, P), P);
  size_t vertexIndex2 = CGAL::test::face_vertex_index(currentFace2, rootSearchVertex2, P);
  Barycentric_coordinate baryCoord2(vertexIndex2 == 0 ? FT(1.0) : FT(0.0), vertexIndex2 == 1 ? FT(1.0) : FT(0.0), vertexIndex2 == 2 ? FT(1.0) : FT(0.0));
  
  std::vector<Polyhedron_shortest_path::Face_location> faceLocations;
  faceLocations.push_back(Polyhedron_shortest_path::Face_location(currentFace, baryCoord));
  faceLocations.push_back(Polyhedron_shortest_path::Face_location(currentFace2, baryCoord2));
  
  shortestPaths.compute_shortest_paths(faceLocations.begin(), faceLocations.end());
  
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
  
  currentVertex = startVertex;
  
  for (size_t i = 0; i < 8; ++i)
  {
    BOOST_CHECK_CLOSE(shortestPaths.shortest_distance_to_vertex(*currentVertex).first, expectedDistances2[i], Kernel::FT(0.0001));
    ++currentVertex;
  }
}
  
BOOST_AUTO_TEST_CASE( test_boundary_mesh )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
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

  Traits traits;
  
  Traits::Project_triangle_3_to_triangle_2 project_triangle_3_to_triangle_2(traits.project_triangle_3_to_triangle_2_object());
  Traits::Compute_squared_distance_3 compute_squared_distance_3(traits.compute_squared_distance_3_object());
  Traits::Compute_squared_distance_2 compute_squared_distance_2(traits.compute_squared_distance_2_object());
  Traits::Construct_barycenter_3 construct_barycenter_3(traits.construct_barycenter_3_object());
  Traits::Flatten_triangle_3_along_segment_2 flatten_triangle_3_along_segment_2(traits.flatten_triangle_3_along_segment_2_object());
  
  struct Construct_barycenter_in_triangle_3
  {
    Traits::Construct_barycenter_3 m_cb3;
  
    Construct_barycenter_in_triangle_3(Traits::Construct_barycenter_3 cb3)
      : m_cb3(cb3)
    {
    }
  
    Point_3 operator() (const Triangle_3& t, const Barycentric_coordinate& b)
    {
      return m_cb3(t[0], b[0], t[1], b[1], t[2], b[2]);
    }
  } construct_barycenter_in_triangle_3(construct_barycenter_3);
  
  std::ifstream inFile("data/boundary_mesh.off");
  
  Polyhedron_3 P;
  
  inFile >> P;
  
  inFile.close();
  
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
  
  Barycentric_coordinate startLocation(FT(0.1), FT(0.8), FT(0.1));
  
  typedef boost::property_map<Polyhedron_3, CGAL::face_external_index_t>::type FaceIndexMap;
  
  FaceIndexMap faceIndexMap(CGAL::get(CGAL::face_external_index, P));
  
  face_descriptor face = *startFace;

  Polyhedron_shortest_path shortestPaths(P, traits);
  //shortestPaths.m_debugOutput = true;
  shortestPaths.compute_shortest_paths(*startFace, startLocation);
  
  Triangle_3 firstTriangle(vertexLocations[1], vertexLocations[0], vertexLocations[2]);
  
  Point_3 locationInTriangle(construct_barycenter_in_triangle_3(firstTriangle, startLocation));
  
  FT dist0 = shortestPaths.shortest_distance_to_vertex(vertexHandles[0]).first;
  BOOST_CHECK_CLOSE(dist0, CGAL::sqrt(compute_squared_distance_3(locationInTriangle, vertexLocations[0])), FT(0.000001));
  
  FT dist1 = shortestPaths.shortest_distance_to_vertex(vertexHandles[1]).first;
  BOOST_CHECK_CLOSE(dist1, CGAL::sqrt(compute_squared_distance_3(locationInTriangle, vertexLocations[1])), FT(0.000001));
  
  FT dist2 = shortestPaths.shortest_distance_to_vertex(vertexHandles[2]).first;
  BOOST_CHECK_CLOSE(dist2, CGAL::sqrt(compute_squared_distance_3(locationInTriangle, vertexLocations[2])), FT(0.000001));
  
  FT dist3 = shortestPaths.shortest_distance_to_vertex(vertexHandles[3]).first;
  BOOST_CHECK_CLOSE(dist3, CGAL::sqrt(compute_squared_distance_3(locationInTriangle, vertexLocations[3])), FT(0.000001));
  
  FT dist4 = shortestPaths.shortest_distance_to_vertex(vertexHandles[4]).first;
  BOOST_CHECK_CLOSE(dist4, CGAL::sqrt(compute_squared_distance_3(locationInTriangle, vertexLocations[1])) + CGAL::sqrt(compute_squared_distance_3(vertexLocations[1], vertexLocations[4])), FT(0.000001));

  FT dist5 = shortestPaths.shortest_distance_to_vertex(vertexHandles[5]).first;
  BOOST_CHECK_CLOSE(dist5, CGAL::sqrt(compute_squared_distance_3(locationInTriangle, vertexLocations[3])) + CGAL::sqrt(compute_squared_distance_3(vertexLocations[3], vertexLocations[5])), FT(0.000001));

  Barycentric_coordinate somewhereElseInFirstTriangle(0.8, 0.05, 0.15);
  
  FT distT0 = shortestPaths.shortest_distance_to_location(faceHandles[0], somewhereElseInFirstTriangle).first;
  BOOST_CHECK_CLOSE(distT0, CGAL::sqrt(compute_squared_distance_3(locationInTriangle, construct_barycenter_in_triangle_3(firstTriangle, somewhereElseInFirstTriangle))), FT(0.000001));
  
  Triangle_3 oneStepTriangle(vertexLocations[4], vertexLocations[1], vertexLocations[3]);
  Barycentric_coordinate locationInOneStepTriangle(0.1, 0.8, 0.1);
  
  CGAL::test::Edge_sequence_collector<Traits> collector(P);
  shortestPaths.shortest_path_sequence(faceHandles[2], locationInOneStepTriangle, collector);

  FT distT2 = shortestPaths.shortest_distance_to_location(faceHandles[2], locationInOneStepTriangle).first;
  BOOST_CHECK_CLOSE(distT2, dist1 + CGAL::sqrt(compute_squared_distance_3(vertexLocations[1], construct_barycenter_in_triangle_3(oneStepTriangle, locationInOneStepTriangle))), FT(0.00001));

  Triangle_3 twoStepTriangle(vertexLocations[6], vertexLocations[5], vertexLocations[7]);
  Barycentric_coordinate locationInTwoStepTriangle(0.8, 0.1, 0.1);
  
  FT distT5 = shortestPaths.shortest_distance_to_location(faceHandles[5], locationInTwoStepTriangle).first;
  BOOST_CHECK_CLOSE(distT5, dist3 + CGAL::sqrt(compute_squared_distance_3(vertexLocations[3], construct_barycenter_in_triangle_3(twoStepTriangle, locationInTwoStepTriangle))), FT(0.00001));
  
  Triangle_3 threeStepTriangle(vertexLocations[7], vertexLocations[5], vertexLocations[8]);
  Barycentric_coordinate locationInThreeStepTriangle(0.2, 0.6, 0.2);
  
  FT distT6 = shortestPaths.shortest_distance_to_location(faceHandles[6], locationInThreeStepTriangle).first;
  BOOST_CHECK_CLOSE(distT6, dist5 + CGAL::sqrt(compute_squared_distance_3(vertexLocations[5], construct_barycenter_in_triangle_3(threeStepTriangle, locationInThreeStepTriangle))), FT(0.00001));
}


// Hack to trick cgal_create_CMakeLists into using this file even without a main
// int main(int argc, char** argv)
// {
//   return 0;
// }
