// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path_traits.h>
#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>

#include <CGAL/Random.h>

#include <CGAL/test_util.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <cmath>

#define BOOST_TEST_MODULE polyhedron_shortest_path_test_3
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE( test_find_nearest_face_location_on_surface )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;
  typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Polyhedron_3> Traits;
  typedef Traits::Barycentric_coordinate Barycentric_coordinate;
  typedef Traits::FT FT;
  typedef Traits::Point_3 Point_3;
  typedef Traits::Triangle_3 Triangle_3;
  typedef boost::graph_traits<Polyhedron_3> GraphTraits;
  typedef GraphTraits::face_descriptor face_descriptor;
  typedef GraphTraits::face_iterator face_iterator;
  typedef CGAL::Polyhedron_shortest_path<Traits> Polyhedron_shortest_path;
  typedef boost::property_map<Polyhedron_3, CGAL::vertex_point_t>::type VPM;
  typedef boost::property_map<Polyhedron_3, CGAL::face_external_index_t>::type FIM;
  
  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron_3, VPM> AABB_face_graph_primitive;
  typedef CGAL::AABB_traits<Kernel, AABB_face_graph_primitive> AABB_face_graph_traits;
  typedef CGAL::AABB_tree<AABB_face_graph_traits> AABB_face_graph_tree;
  
  Traits traits;
  
  Traits::Construct_barycenter_3 construct_barycenter_3(traits.construct_barycenter_3_object());
  
  Polyhedron_3 polyhedron;
  
  std::ifstream in("data/anchor.off");
    
  in >> polyhedron;
  
  in.close();

  Polyhedron_shortest_path shortestPaths(polyhedron, traits);
  
  face_iterator facesBegin, facesEnd;
  boost::tie(facesBegin, facesEnd) = CGAL::faces(polyhedron);
  
  std::vector<face_descriptor> facesList;
  
  for (face_iterator facesCurrent = facesBegin; facesCurrent != facesEnd; ++facesCurrent)
  {
    facesList.push_back(*facesCurrent);
  }
  
  CGAL::Random random(6008991);
  
  size_t numTrials = 30;

  FIM faceIndexMap(CGAL::get(CGAL::face_external_index, polyhedron));
  VPM vertexPointMap(CGAL::get(CGAL::vertex_point, polyhedron));
 
  for (size_t i = 0; i < numTrials; ++i)
  {
    size_t faceIndex = random.get_int(0, facesList.size());
    face_descriptor face = facesList[faceIndex];

    Triangle_3 faceTriangle = CGAL::internal::triangle_from_halfedge<Triangle_3, Polyhedron_3, VPM>(CGAL::halfedge(face, polyhedron), polyhedron, vertexPointMap);
    
    Barycentric_coordinate location = CGAL::test::random_coordinate<Traits>(random);
    
    Point_3 location3d = construct_barycenter_3(faceTriangle[0], location[0], faceTriangle[1], location[1], faceTriangle[2], location[2]);
    
    Polyhedron_shortest_path::Face_location faceLocation = shortestPaths.locate<AABB_face_graph_traits>(location3d);
    
    BOOST_CHECK_EQUAL(faceIndexMap[face], faceIndexMap[faceLocation.first]);
    BOOST_CHECK_CLOSE(location[0], faceLocation.second[0], FT(0.0001));
    BOOST_CHECK_CLOSE(location[1], faceLocation.second[1], FT(0.0001));
    BOOST_CHECK_CLOSE(location[2], faceLocation.second[2], FT(0.0001));
  }
}

// Hack to trick cgal_create_CMakeLists into using this file even without a main
// int main(int argc, char** argv)
// {
//   return 0;
// }
