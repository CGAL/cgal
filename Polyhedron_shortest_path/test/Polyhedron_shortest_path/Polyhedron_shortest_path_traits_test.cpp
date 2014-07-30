// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path_traits.h>
#include <CGAL/Polyhedron_shortest_path/internal/function_objects.h>
#include <CGAL/Polyhedron_shortest_path/internal/Barycentric.h>
#include <CGAL/Polyhedron_shortest_path/internal/misc_functions.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/test_util.h>

#include <iostream>
#include <fstream>
#include <utility>

#define BOOST_TEST_MODULE polyhedron_shortest_path_traits_test
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE( project_triangle3D_to_triangle2D )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
  typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Polyhedron_3> Traits;
  
  Traits traits;
  Traits::Compute_squared_distance_2 compute_squared_distance_2(traits.compute_squared_distance_2_object());
  Traits::Project_triangle_3_to_triangle_2 project_triangle_3_to_triangle_2(traits.project_triangle_3_to_triangle_2_object());
  
  Traits::Triangle_3 sourceTriangle(
    Traits::Point_3(Kernel::FT(0), Kernel::FT(1), Kernel::FT(2)), 
    Traits::Point_3(Kernel::FT(5), Kernel::FT(4), Kernel::FT(3)),
    Traits::Point_3(Kernel::FT(4), Kernel::FT(9), Kernel::FT(16)));

  Traits::Triangle_2 transformedTriangle = project_triangle_3_to_triangle_2(sourceTriangle);
  BOOST_CHECK_CLOSE(compute_squared_distance_2(sourceTriangle[1], sourceTriangle[0]), compute_squared_distance_2(transformedTriangle[1], transformedTriangle[0]), Kernel::FT(0.0000001));
  BOOST_CHECK_CLOSE(compute_squared_distance_2(sourceTriangle[2], sourceTriangle[0]), compute_squared_distance_2(transformedTriangle[2], transformedTriangle[0]), Kernel::FT(0.0000001));
  BOOST_CHECK_CLOSE(compute_squared_distance_2(sourceTriangle[2], sourceTriangle[1]), compute_squared_distance_2(transformedTriangle[2], transformedTriangle[1]), Kernel::FT(0.0000001));
  
}

BOOST_AUTO_TEST_CASE( test_simple_2D_barycentric_coordinates )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
  typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Polyhedron_3> Traits;

  Traits traits;
  
  Traits::Construct_barycentric_coordinate_in_triangle_2 construct_barycentric_coordinate_in_triangle_2(traits.construct_barycentric_coordinate_in_triangle_2_object());
  Traits::Construct_barycenter_2 construct_barycenter_2(traits.construct_barycenter_2_object());
  Traits::Classify_barycentric_coordinate classify_barycentric_coordinate(traits.classify_barycentric_coordinate_object());
  
  // effectively a 1-1 mapping triangle for barycentric coords
  Traits::Triangle_2 simpleTriangle(
    Traits::Point_2(Kernel::FT(0), Kernel::FT(0)), 
    Traits::Point_2(Kernel::FT(1), Kernel::FT(0)),
    Traits::Point_2(Kernel::FT(0), Kernel::FT(1)));

  Traits::Barycentric_coordinate b0 = construct_barycentric_coordinate_in_triangle_2(simpleTriangle, simpleTriangle[0]);
  
  BOOST_CHECK_CLOSE(Kernel::FT(1.0), b0[0], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(0.0), b0[1], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(0.0), b0[2], Kernel::FT(0.000001));
  
  size_t outVertex0;
  CGAL::internal::Barycentric_coordinate_type b0Type;
  boost::tie(b0Type, outVertex0) = classify_barycentric_coordinate(b0);
  
  BOOST_CHECK_EQUAL(b0Type, CGAL::internal::BARYCENTRIC_COORDINATE_VERTEX);
  BOOST_CHECK_EQUAL(outVertex0, 0);

  Traits::Barycentric_coordinate b1 = construct_barycentric_coordinate_in_triangle_2(simpleTriangle, simpleTriangle[1]);
  
  BOOST_CHECK_CLOSE(Kernel::FT(0.0), b1[0], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(1.0), b1[1], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(0.0), b1[2], Kernel::FT(0.000001));
  
  size_t outVertex1;
  CGAL::internal::Barycentric_coordinate_type b1Type;
  boost::tie(b1Type, outVertex1) = classify_barycentric_coordinate(b1);
  
  BOOST_CHECK_EQUAL(b1Type, CGAL::internal::BARYCENTRIC_COORDINATE_VERTEX);
  BOOST_CHECK_EQUAL(outVertex1, 1);
  
  Traits::Barycentric_coordinate b2 = construct_barycentric_coordinate_in_triangle_2(simpleTriangle, simpleTriangle[2]);
  
  BOOST_CHECK_CLOSE(Kernel::FT(0.0), b2[0], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(0.0), b2[1], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(1.0), b2[2], Kernel::FT(0.000001));
  
  size_t outVertex2;
  CGAL::internal::Barycentric_coordinate_type b2Type;
  boost::tie(b2Type, outVertex2) = classify_barycentric_coordinate(b2);
  
  BOOST_CHECK_EQUAL(b2Type, CGAL::internal::BARYCENTRIC_COORDINATE_VERTEX);
  BOOST_CHECK_EQUAL(outVertex2, 2);
  
  Traits::Point_2 location(Kernel::FT(0.3), Kernel::FT(0.2));
  Traits::Barycentric_coordinate bLocation = construct_barycentric_coordinate_in_triangle_2(simpleTriangle, location);
  
  size_t dummyOut;
  CGAL::internal::Barycentric_coordinate_type bLocationType;
  boost::tie(bLocationType, dummyOut) = classify_barycentric_coordinate(bLocation);
  
  BOOST_CHECK_EQUAL(bLocationType, CGAL::internal::BARYCENTRIC_COORDINATE_INTERNAL);
  
  BOOST_CHECK_CLOSE(Kernel::FT(1.0) - location[0] - location[1], bLocation[0], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(location[0], bLocation[1], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(location[1], bLocation[2], Kernel::FT(0.000001));
  
  Traits::Point_2 outLocation = construct_barycenter_2(simpleTriangle[0], bLocation[0], simpleTriangle[1], bLocation[1], simpleTriangle[2], bLocation[2]);
  BOOST_CHECK_CLOSE(location[0], outLocation[0], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(location[1], outLocation[1], Kernel::FT(0.000001));
}
  
BOOST_AUTO_TEST_CASE( barycentric_coords_3D )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
  typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Polyhedron_3> Traits;
  
  Traits traits;
  
  Traits::Construct_barycentric_coordinate_in_triangle_3 construct_barycentric_coordinate_in_triangle_3(traits.construct_barycentric_coordinate_in_triangle_3_object());
  Traits::Construct_barycenter_3 construct_barycenter_3(traits.construct_barycenter_3_object());

  Traits::Triangle_3 quadrantTriangle(
    Traits::Point_3(Kernel::FT(1.0), Kernel::FT(0.0), Kernel::FT(1.0)),
    Traits::Point_3(Kernel::FT(-1.0), Kernel::FT(-1.0), Kernel::FT(-1.0)),
    Traits::Point_3(Kernel::FT(-1.0), Kernel::FT(1.0), Kernel::FT(-1.0)));
    
  Traits::Barycentric_coordinate b0 = construct_barycentric_coordinate_in_triangle_3(quadrantTriangle, quadrantTriangle[0]);
  
  BOOST_CHECK_CLOSE(Kernel::FT(1.0), b0[0], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(0.0), b0[1], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(0.0), b0[2], Kernel::FT(0.000001));
  
  Traits::Barycentric_coordinate b1 = construct_barycentric_coordinate_in_triangle_3(quadrantTriangle, quadrantTriangle[1]);
  
  BOOST_CHECK_CLOSE(Kernel::FT(0.0), b1[0], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(1.0), b1[1], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(0.0), b1[2], Kernel::FT(0.000001));
  
  Traits::Barycentric_coordinate b2 = construct_barycentric_coordinate_in_triangle_3(quadrantTriangle, quadrantTriangle[2]);
  
  BOOST_CHECK_CLOSE(Kernel::FT(0.0), b2[0], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(0.0), b2[1], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(1.0), b2[2], Kernel::FT(0.000001));
  
  Traits::Barycentric_coordinate bOrigin = construct_barycentric_coordinate_in_triangle_3(quadrantTriangle, Traits::Point_3(CGAL::ORIGIN));

  BOOST_CHECK_CLOSE(Kernel::FT(0.5), bOrigin[0], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(0.25), bOrigin[1], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(0.25), bOrigin[2], Kernel::FT(0.000001));
  
  Traits::Point_3 originOutAgain = construct_barycenter_3(quadrantTriangle[0], bOrigin[0], quadrantTriangle[1], bOrigin[1], quadrantTriangle[2], bOrigin[2]);
  
  BOOST_CHECK_CLOSE(Kernel::FT(0.0), originOutAgain[0], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(0.0), originOutAgain[1], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(0.0), originOutAgain[2], Kernel::FT(0.000001));
  
  Traits::Point_3 pNegative(Kernel::FT(-0.5), Kernel::FT(-0.5), Kernel::FT(-0.5));
  Traits::Barycentric_coordinate bNegative = construct_barycentric_coordinate_in_triangle_3(quadrantTriangle, pNegative);

  BOOST_CHECK_CLOSE(Kernel::FT(0.25), bNegative[0], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(0.625), bNegative[1], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(0.125), bNegative[2], Kernel::FT(0.000001));
  
  Traits::Point_3 negativeOutAgain = construct_barycenter_3(quadrantTriangle[0], bNegative[0], quadrantTriangle[1], bNegative[1], quadrantTriangle[2], bNegative[2]);
  
  BOOST_CHECK_CLOSE(Kernel::FT(-0.5), negativeOutAgain[0], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(-0.5), negativeOutAgain[1], Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(-0.5), negativeOutAgain[2], Kernel::FT(0.000001));
}

BOOST_AUTO_TEST_CASE( simple_flattening_triangle_along_edge )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
  typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Polyhedron_3> Traits;
  
  Traits traits;
  
  Traits::Flatten_triangle_3_along_segment_2 flatten_triangle_3_along_segment_2(traits.flatten_triangle_3_along_segment_2_object());
  Traits::Compute_squared_distance_2 compute_squared_distance_2(traits.compute_squared_distance_2_object());
  Traits::Compute_squared_distance_3 compute_squared_distance_3(traits.compute_squared_distance_3_object());
  Traits::Orientation_2 orientation_2(traits.orientation_2_object());
  
  const Kernel::FT baseDistance(10.0);
  const Kernel::FT height(6.0);

  Traits::Triangle_3 sourceTriangle(
    Traits::Point_3(Kernel::FT(-4), Kernel::FT(0), Kernel::FT(3)), 
    Traits::Point_3(Kernel::FT(4), Kernel::FT(0), Kernel::FT(-3)),
    Traits::Point_3(Kernel::FT(0), Kernel::FT(height), Kernel::FT(0)));
    
  Traits::Segment_2 layoutEdge(Traits::Point_2(CGAL::ORIGIN), Traits::Point_2(baseDistance, Kernel::FT(0.0)));
  
  Traits::Triangle_2 flattened = flatten_triangle_3_along_segment_2(sourceTriangle, 0, layoutEdge);

  BOOST_CHECK_CLOSE(compute_squared_distance_3(sourceTriangle[0], sourceTriangle[1]), Kernel::FT(100.0), Kernel::FT(0.000001));
  BOOST_CHECK_EQUAL(orientation_2(flattened[0], flattened[1], flattened[2]), CGAL::LEFT_TURN);
  
  BOOST_CHECK_CLOSE(Kernel::FT(flattened[0].x()), layoutEdge[0].x(), Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(flattened[0].y()), layoutEdge[0].y(), Kernel::FT(0.000001));
  
  BOOST_CHECK_CLOSE(Kernel::FT(flattened[1].x()), layoutEdge[1].x(), Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(flattened[1].y()), layoutEdge[1].y(), Kernel::FT(0.000001));

  BOOST_CHECK_CLOSE(Kernel::FT(flattened[2].x()), baseDistance * (-sourceTriangle[0].x() / (sourceTriangle[1].x() - sourceTriangle[0].x())), Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(Kernel::FT(flattened[2].y()), height, Kernel::FT(0.000001));
  
  BOOST_CHECK_CLOSE(compute_squared_distance_3(sourceTriangle[0], sourceTriangle[1]), compute_squared_distance_2(flattened[0], flattened[1]), Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(compute_squared_distance_3(sourceTriangle[0], sourceTriangle[2]), compute_squared_distance_2(flattened[0], flattened[2]), Kernel::FT(0.000001));
  BOOST_CHECK_CLOSE(compute_squared_distance_3(sourceTriangle[1], sourceTriangle[2]), compute_squared_distance_2(flattened[1], flattened[2]), Kernel::FT(0.000001));
}

BOOST_AUTO_TEST_CASE( nonsimple_flattening_triangle_along_edge )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
  typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Polyhedron_3> Traits;
  
  Traits traits;
  
  Traits::Flatten_triangle_3_along_segment_2 flatten_triangle_3_along_segment_2(traits.flatten_triangle_3_along_segment_2_object());
  Traits::Compute_squared_distance_2 compute_squared_distance_2(traits.compute_squared_distance_2_object());
  Traits::Compute_squared_distance_3 compute_squared_distance_3(traits.compute_squared_distance_3_object());
  Traits::Orientation_2 orientation_2(traits.orientation_2_object());

  Traits::Triangle_3 sourceTriangle(
    Traits::Point_3(Kernel::FT(-5), Kernel::FT(16), Kernel::FT(-3)), 
    Traits::Point_3(Kernel::FT(5), Kernel::FT(-9), Kernel::FT(7)),
    Traits::Point_3(Kernel::FT(0), Kernel::FT(4), Kernel::FT(5)));
    
  for (size_t edgeIndex = 0; edgeIndex < 3; ++edgeIndex)
  {    
    const Kernel::FT baseDistance = CGAL::sqrt(compute_squared_distance_3(sourceTriangle.vertex(edgeIndex), sourceTriangle.vertex(edgeIndex + 1)));
    const Traits::Vector_2 direction(Kernel::FT(3.0) / Kernel::FT(5.0), Kernel::FT(4.0) / Kernel::FT(5.0));  
    
    Traits::Segment_2 layoutEdge(Traits::Point_2(CGAL::ORIGIN), Traits::Point_2(CGAL::ORIGIN) + direction * baseDistance);

    Traits::Triangle_2 flattened = flatten_triangle_3_along_segment_2(sourceTriangle, edgeIndex, layoutEdge);
    
    BOOST_CHECK_EQUAL(orientation_2(flattened.vertex(edgeIndex), flattened.vertex(edgeIndex + 1), flattened.vertex(edgeIndex + 2)), CGAL::LEFT_TURN);
    
    BOOST_CHECK_CLOSE(Kernel::FT(flattened.vertex(edgeIndex).x()), layoutEdge[0].x(), Kernel::FT(0.000001));
    BOOST_CHECK_CLOSE(Kernel::FT(flattened.vertex(edgeIndex).y()), layoutEdge[0].y(), Kernel::FT(0.000001));
    
    BOOST_CHECK_CLOSE(Kernel::FT(flattened.vertex(edgeIndex + 1).x()), layoutEdge[1].x(), Kernel::FT(0.000001));
    BOOST_CHECK_CLOSE(Kernel::FT(flattened.vertex(edgeIndex + 1).y()), layoutEdge[1].y(), Kernel::FT(0.000001));
    
    BOOST_CHECK_CLOSE(compute_squared_distance_3(sourceTriangle[0], sourceTriangle[1]), compute_squared_distance_2(flattened[0], flattened[1]), Kernel::FT(0.000001));
    BOOST_CHECK_CLOSE(compute_squared_distance_3(sourceTriangle[0], sourceTriangle[2]), compute_squared_distance_2(flattened[0], flattened[2]), Kernel::FT(0.000001));
    BOOST_CHECK_CLOSE(compute_squared_distance_3(sourceTriangle[1], sourceTriangle[2]), compute_squared_distance_2(flattened[1], flattened[2]), Kernel::FT(0.000001));
  }
}
  
BOOST_AUTO_TEST_CASE( detect_is_saddle_vertex )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
  typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Polyhedron_3> Traits;
  
  typedef boost::graph_traits<Polyhedron_3> GraphTraits;
  typedef GraphTraits::vertex_descriptor vertex_descriptor;
  typedef GraphTraits::vertex_iterator vertex_iterator;
  
  Traits traits;
  Traits::Is_saddle_vertex is_saddle_vertex(traits.is_saddle_vertex_object());

  std::ifstream inFile("data/saddle_vertex_mesh.off");
  
  Polyhedron_3 P;
  
  inFile >> P;
  
  inFile.close();
  
  size_t currentVertexIndex = 0;
  
  vertex_iterator currentVertex;
  vertex_iterator endVertex;

  for (boost::tie(currentVertex, endVertex) = boost::vertices(P); currentVertex != endVertex; ++currentVertex)
  {
    if (currentVertexIndex <= 3 || currentVertexIndex == 7)
    {
      BOOST_CHECK(!is_saddle_vertex(*currentVertex, P));
    }
    else
    {
      BOOST_CHECK(is_saddle_vertex(*currentVertex, P));
    }

    ++currentVertexIndex;
  }
}

// Hack to trick cgal_test_with_cmake into using this file even without a main
// int main(int argc, char** argv)
// {
//   return 0;
// }
