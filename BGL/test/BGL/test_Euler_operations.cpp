#define BOOST_TEST_MODULE euler_operations test
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>

#include "test_Prefix.h"
#include <CGAL/boost/graph/Euler_operations.h>
#include <boost/range/algorithm.hpp>


BOOST_AUTO_TEST_CASE_TEMPLATE( join_face_test, T, test_graphs )
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);

  Surface_fixture_1<T> f;

  bool found;
  halfedge_descriptor e;
  boost::tie(e, found) = halfedge(f.w, f.v, f.m);
  BOOST_CHECK(found);
  // manually set the halfedge of f.f1 to the edge that is to be
  // removed to provoke a special case
  set_halfedge(f.f1, e, f.m);
  CGAL::Euler::join_face(e,f.m);

  BOOST_CHECK(CGAL::internal::exact_num_faces(f.m) == 2);
  BOOST_CHECK(CGAL::internal::exact_num_edges(f.m) == 6);
  
  CGAL::Halfedge_around_face_iterator<T> begin, end;
  boost::tie(begin, end) = halfedges_around_face(halfedge(f.f1, f.m), f.m);
  BOOST_CHECK(std::distance(begin, end) == 4);
  for(; begin != end; ++begin)
  {

    halfedge_descriptor hd = *begin;
    BOOST_CHECK(face(hd, f.m) == f.f1);

  }

  face_iterator fit, fend;
  for(boost::tie(fit, fend) = faces(f.m); fit != fend; ++fit) {
    BOOST_CHECK(*fit == f.f1 || *fit == f.f3);
  }
  
  BOOST_CHECK(degree(f.w, f.m) == 2);
  BOOST_CHECK(degree(f.v, f.m) == 3);
  BOOST_CHECK(CGAL::is_valid(f.m));

}



BOOST_AUTO_TEST_CASE_TEMPLATE( remove_face_test_1, T, test_graphs )
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);
  Surface_fixture_1<T> f;

  // find the edge between x and y
  bool found;
  halfedge_descriptor e;
  boost::tie(e, found) = halfedge(f.x, f.y, f.m);
  BOOST_CHECK(found);
  BOOST_CHECK(face(e, f.m) == f.f3);

  CGAL::Euler::remove_face(e,f.m);

  BOOST_CHECK(CGAL::is_valid(f.m));

  BOOST_CHECK_EQUAL(degree(f.v, f.m), 3);
  BOOST_CHECK_EQUAL(degree(f.x, f.m), 2);
  BOOST_CHECK_EQUAL(CGAL::internal::exact_num_faces(f.m), 2);
  BOOST_CHECK_EQUAL(CGAL::internal::exact_num_edges(f.m), 5);
  BOOST_CHECK_EQUAL(CGAL::internal::exact_num_vertices(f.m), 4);
  halfedge_iterator eb, ee;
  int count = 0;
  for(boost::tie(eb, ee) = halfedges(f.m); eb != ee; ++eb) {
    if(face(*eb,f.m) == boost::graph_traits<T>::null_face())
      ++count;
  }
  BOOST_CHECK_EQUAL(count, 4);
}



BOOST_AUTO_TEST_CASE_TEMPLATE( remove_face_test_2, T, test_graphs )
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);

  Surface_fixture_2<T> f;
 
  // find the edge between x and v
  bool found;
  halfedge_descriptor e;

  boost::tie(e, found) = halfedge(f.x, f.w, f.m);
  BOOST_CHECK(found);
  boost::tie(e, found) = halfedge(f.x, f.v, f.m);
  BOOST_CHECK(found);
  BOOST_CHECK(face(e, f.m) == f.f1);
  CGAL::Euler::remove_face(e,f.m);
  BOOST_CHECK(CGAL::is_valid(f.m));

  BOOST_CHECK(CGAL::internal::exact_num_faces(f.m) == 3);
  BOOST_CHECK(CGAL::internal::exact_num_edges(f.m) == 7);
  BOOST_CHECK(CGAL::internal::exact_num_vertices(f.m) == 5);

  boost::tie(e, found) = halfedge(f.x, f.w, f.m);
  BOOST_CHECK(found);
  BOOST_CHECK(face(e,f.m) == boost::graph_traits<T>::null_face());
  
  // check the boundary

  halfedge_descriptor n = next(e, f.m);
  while(n != e) {
    BOOST_CHECK(face(n,f.m) == boost::graph_traits<T>::null_face() );
    n = next(n, f.m); 
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE( add_face_to_border_test, T, test_graphs )
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);
 
  Surface_fixture_5<T> f;

  CGAL::Euler::add_face_to_border(f.h1, f.h2, f.m);

  BOOST_CHECK(CGAL::is_valid(f.m));

}


BOOST_AUTO_TEST_CASE_TEMPLATE( join_vertex_interior_test, T, test_graphs )
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);

  Surface_fixture_3<T> f;
  halfedge_descriptor e;

  bool found;
  boost::tie(e, found) = halfedge(f.w, f.x, f.m);
  BOOST_CHECK(found);
  CGAL::Euler::join_vertex(e,f.m);
  BOOST_CHECK_EQUAL(CGAL::internal::exact_num_faces(f.m), 2);
  BOOST_CHECK_EQUAL(CGAL::internal::exact_num_vertices(f.m), 5);
  BOOST_CHECK_EQUAL(CGAL::internal::exact_num_edges(f.m), 6);
  BOOST_CHECK_EQUAL(boost::distance(halfedges_around_face(halfedge(f.f1, f.m), f.m)), 3);
  BOOST_CHECK_EQUAL(boost::distance(halfedges_around_face(halfedge(f.f2, f.m), f.m)), 3);
  BOOST_CHECK_EQUAL(degree(f.x, f.m), 4);
  BOOST_CHECK(CGAL::is_valid(f.m));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( join_vertex_exterior_test, T, test_graphs )
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);

  {
    // exterior edge is collapsed
    Surface_fixture_3<T> f;
    halfedge_descriptor e;
    bool found;
    boost::tie(e, found) = halfedge(f.w, f.y, f.m);
    assert(source(e,f.m) == f.w);
    assert(target(e,f.m) == f.y);
    BOOST_CHECK(found);
    CGAL::Euler::join_vertex(e,f.m);
    BOOST_CHECK_EQUAL(CGAL::internal::exact_num_faces(f.m), 2);
    BOOST_CHECK_EQUAL(CGAL::internal::exact_num_vertices(f.m), 5);
    BOOST_CHECK_EQUAL(CGAL::internal::exact_num_edges(f.m), 6);
    BOOST_CHECK_EQUAL(boost::distance(halfedges_around_face(halfedge(f.f1, f.m), f.m)), 4);
    BOOST_CHECK_EQUAL(boost::distance(halfedges_around_face(halfedge(f.f2, f.m), f.m)), 3);
    BOOST_CHECK_EQUAL(degree(f.y, f.m), 3);
    BOOST_CHECK(CGAL::is_valid(f.m));
  }

  {
    Surface_fixture_3<T> f;
    halfedge_descriptor e;
    bool found;
    boost::tie(e, found) = halfedge(f.y, f.w, f.m);

    assert(source(e,f.m) == f.y);
    assert(target(e,f.m) == f.w);
    BOOST_CHECK(found);
    CGAL::Euler::join_vertex(e,f.m);
    BOOST_CHECK_EQUAL(CGAL::internal::exact_num_faces(f.m), 2);
    BOOST_CHECK_EQUAL(CGAL::internal::exact_num_vertices(f.m), 5);
    BOOST_CHECK_EQUAL(CGAL::internal::exact_num_edges(f.m), 6);
    BOOST_CHECK_EQUAL(boost::distance(halfedges_around_face(halfedge(f.f1, f.m), f.m)), 4);
    BOOST_CHECK_EQUAL(boost::distance(halfedges_around_face(halfedge(f.f2, f.m), f.m)), 3);
 
    BOOST_CHECK(CGAL::is_valid(f.m));
    BOOST_CHECK_EQUAL(degree(f.w, f.m), 3);

  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE( split_vertex, T, test_graphs )
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);

  Surface_fixture_3<T> f;
  halfedge_descriptor h1, h2;
  bool found;
  boost::tie(h1, found) = halfedge(f.w, f.y, f.m);
  BOOST_CHECK(found);
  boost::tie(h2, found) = halfedge(f.z, f.y, f.m);
  BOOST_CHECK(found);
  BOOST_CHECK(face(h2, f.m) == Traits::null_face());

  // split border vertex y
  CGAL::Euler::split_vertex(h1, h2,f.m);
  BOOST_CHECK(CGAL::is_valid(f.m));
  BOOST_CHECK_EQUAL(CGAL::internal::exact_num_vertices(f.m), 7);
  BOOST_CHECK_EQUAL(CGAL::internal::exact_num_edges(f.m), 8);
  BOOST_CHECK_EQUAL(boost::distance(halfedges_around_face(h1, f.m)), 5);
  BOOST_CHECK_EQUAL(boost::distance(halfedges_around_face(h2, f.m)), 7);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( split_join_vertex_inverse, T, test_graphs )
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);
  Surface_fixture_3<T> f;
  halfedge_descriptor h, h1, h2;
  bool found;
  boost::tie(h, found) = halfedge(f.w, f.x, f.m);
  BOOST_CHECK(found);
  CGAL::Euler::join_vertex(h,f.m);
  BOOST_CHECK(CGAL::is_valid(f.m));
  boost::tie(h1, found) = halfedge(f.z, f.x, f.m);
  BOOST_CHECK(found);
  boost::tie(h2, found) = halfedge(f.v, f.x, f.m);
  BOOST_CHECK(found);
  CGAL::Euler::join_vertex(CGAL::Euler::split_vertex(h1, h2,f.m),f.m);
  BOOST_CHECK(CGAL::is_valid(f.m));
  
  BOOST_CHECK_EQUAL(CGAL::internal::exact_num_vertices(f.m), 5);
  BOOST_CHECK_EQUAL(CGAL::internal::exact_num_faces(f.m), 2);
  BOOST_CHECK_EQUAL(CGAL::internal::exact_num_edges(f.m), 6);
  BOOST_CHECK_EQUAL(CGAL::internal::exact_num_halfedges(f.m), 12);
  BOOST_CHECK_EQUAL(boost::distance(halfedges_around_face(h1, f.m)), 3);
  BOOST_CHECK_EQUAL(boost::distance(halfedges_around_face(h2, f.m)), 3);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( join_loop_test, T, test_graphs )
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);
  Surface_fixture_4<T> f;

  CGAL::Euler::join_loop(f.h1, f.h2, f.m);
  
  BOOST_CHECK(CGAL::is_valid(f.m));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( split_loop_test, T, test_graphs )
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);
  Surface_fixture_8<T> f;

  CGAL::Euler::split_loop(f.h1, f.h2, f.h3, f.m);
  BOOST_CHECK_EQUAL(CGAL::internal::exact_num_vertices(f.m), 8);
  BOOST_CHECK_EQUAL(CGAL::internal::exact_num_faces(f.m), 8);
  BOOST_CHECK_EQUAL(CGAL::internal::exact_num_halfedges(f.m), 24);
  BOOST_CHECK(CGAL::is_valid(f.m));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( split_face_test, T, test_graphs )
{
 Surface_fixture_6<T> f;
 CGAL::Euler::split_face(f.h1, f.h2,f.m);
 BOOST_CHECK_EQUAL(num_vertices(f.m), 4);
 BOOST_CHECK_EQUAL(num_faces(f.m), 2);
 BOOST_CHECK_EQUAL(num_halfedges(f.m), 10);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( make_hole_test, T, test_graphs )
{
 Surface_fixture_7<T> f;
 std::size_t nv = num_vertices(f.m);
 std::size_t nf = num_faces(f.m);
 std::size_t nh = num_halfedges(f.m);

 CGAL::Euler::make_hole(f.h, f.m);

 BOOST_CHECK_EQUAL(CGAL::internal::exact_num_vertices(f.m), nv);
 BOOST_CHECK_EQUAL(CGAL::internal::exact_num_faces(f.m), nf-1 );
 BOOST_CHECK_EQUAL(CGAL::internal::exact_num_halfedges(f.m), nh);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( remove_center_vertex_test, T, test_graphs )
{
 Surface_fixture_7<T> f;
 std::size_t nv = num_vertices(f.m);
 std::size_t nf = num_faces(f.m);
 std::size_t nh = num_halfedges(f.m);

 typename boost::graph_traits<T>::degree_size_type deg = degree(target(f.h,f.m),f.m);
 CGAL::Euler::remove_center_vertex(f.h,f.m);

 BOOST_CHECK_EQUAL(CGAL::internal::exact_num_vertices(f.m), nv-1);
 BOOST_CHECK_EQUAL(CGAL::internal::exact_num_faces(f.m), (nf-deg)+1);
 BOOST_CHECK_EQUAL(CGAL::internal::exact_num_halfedges(f.m), nh-(2*deg));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( join_split_inverse, T, test_graphs )
{
  
}

BOOST_AUTO_TEST_CASE_TEMPLATE( satisfies_link_condition, T, test_graphs )
{
  Surface_fixture_7<T> f;

  BOOST_CHECK(CGAL::Euler::safisfies_link_condition(*edges(f.m).first,f.m));
}


// trick cgal_test_with_cmake
// int main()
// {
//   return 0;
// }
