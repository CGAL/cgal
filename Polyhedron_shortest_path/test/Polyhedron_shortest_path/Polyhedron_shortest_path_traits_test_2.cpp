// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#include <iomanip>

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
#include <cmath>

#define BOOST_TEST_MODULE polyhedron_shortest_path_traits_test_2
#include <boost/test/included/unit_test.hpp>

/*
BOOST_AUTO_TEST_CASE( test_orientation_weirdness )
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
  typedef boost::property_map<typename Traits::Polyhedron, boost::vertex_external_index_t>::type VIM;
  typedef boost::property_map<typename Traits::Polyhedron, boost::edge_external_index_t>::type EIM;
  typedef boost::property_map<typename Traits::Polyhedron, CGAL::halfedge_external_index_t>::type HIM;
  typedef boost::property_map<typename Traits::Polyhedron, CGAL::face_external_index_t>::type FIM;
  
  Traits traits;
  Traits::Orientation_2 orientation_2(traits.orientation_2_object());
  Traits::Intersect_2 intersect_2(traits.intersect_2_object());
  
  std::cout << CGAL::LEFT_TURN << " " << CGAL::RIGHT_TURN << " " << CGAL::COLLINEAR << std::endl;
  
  std::cout << orientation_2(Point_2(FT(0.125437), FT(0.0)), Point_2
  
  orientation_2(sourceImagePoint, m_windowLeft, point) == CGAL::RIGHT_TURN && orientation_2(sourceImagePoint, m_windowRight, point) == CGAL::LEFT_TURN;
}
  */
  
// Hack to trick cgal_test_with_cmake into using this file even without a main
// int main(int argc, char** argv)
// {
//   return 0;
// }