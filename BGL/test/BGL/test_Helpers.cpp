#define BOOST_TEST_MODULE graph_traits_helpers_test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>

#include "test_Prefix.h"

#include <CGAL/BGL/Helper.h>

BOOST_AUTO_TEST_CASE(constant_vertex_is_border_surface_mesh)
{
  Surface_fixture_1<SM> f1;
  Surface_fixture_2<SM> f2;
  Surface_fixture_3<SM> f3;
  BOOST_CHECK(CGAL::internal::constant_vertex_is_border(f1.m));
  BOOST_CHECK(CGAL::internal::constant_vertex_is_border(f2.m));
  BOOST_CHECK(CGAL::internal::constant_vertex_is_border(f3.m));
  
  // invalidate
  boost::graph_traits<SM>::halfedge_descriptor h;
  bool found;
  boost::tie(h, found) = halfedge(f1.w ,f1.v, f1.m);
  BOOST_CHECK(found);

  set_halfedge(f1.v, h, f1.m);
  BOOST_CHECK(!CGAL::internal::constant_vertex_is_border(f1.m));
  CGAL::internal::set_constant_vertex_is_border(f1.m);
  BOOST_CHECK(CGAL::internal::constant_vertex_is_border(f1.m));
}

BOOST_AUTO_TEST_CASE(constant_vertex_is_border_polyhedron)
{
  Surface_fixture_1<Polyhedron> f1;
  Surface_fixture_2<Polyhedron> f2;
  Surface_fixture_3<Polyhedron> f3;
  // those tests are expected to fail for polyhedron, as the file
  // reader does not respect the border guarantee
  BOOST_CHECK(CGAL::internal::constant_vertex_is_border(f1.m));
  BOOST_CHECK(CGAL::internal::constant_vertex_is_border(f2.m));
  BOOST_CHECK(CGAL::internal::constant_vertex_is_border(f3.m));

  // repair them
  CGAL::internal::set_constant_vertex_is_border(f1.m);
  CGAL::internal::set_constant_vertex_is_border(f2.m);
  CGAL::internal::set_constant_vertex_is_border(f3.m);
  
  BOOST_CHECK(f1.m.is_valid());
  BOOST_CHECK(f2.m.is_valid());
  BOOST_CHECK(f3.m.is_valid());

  // check them again
  BOOST_CHECK(CGAL::internal::constant_vertex_is_border(f1.m));
  BOOST_CHECK(CGAL::internal::constant_vertex_is_border(f2.m));
  BOOST_CHECK(CGAL::internal::constant_vertex_is_border(f3.m));
}



// trick cgal_create_CMakeLists
// int main()
// {
//   return 0;
// }
