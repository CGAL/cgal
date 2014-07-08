#define BOOST_TEST_MODULE surface_mesh circulator test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include <CGAL/circulator.h>

#include "SM_common.h"

// trivial circulator testing
BOOST_FIXTURE_TEST_CASE(circ_size, Surface_fixture)
{
  Sm::Vertex_around_target_circulator vvc(m.halfedge(v),m);
  Sm::Halfedge_around_target_circulator hvc(m.halfedge(v),m);
  Sm::Face_around_target_circulator fvc(m.halfedge(v),m);
  Sm::Vertex_around_face_circulator vaf(m.halfedge(f1),m);
  Sm::Halfedge_around_face_circulator hfc(m.halfedge(f1),m);
  
  BOOST_CHECK_EQUAL(CGAL::circulator_size(vvc), 4);
  BOOST_CHECK_EQUAL(CGAL::circulator_size(hvc), 4);
  BOOST_CHECK_EQUAL(CGAL::circulator_size(fvc), 3);
  BOOST_CHECK_EQUAL(CGAL::circulator_size(vaf), 3);
  BOOST_CHECK_EQUAL(CGAL::circulator_size(hfc), 3);
  
}

BOOST_FIXTURE_TEST_CASE(circ_distance, Surface_fixture)
{
  // the distance against end
  Sm::Vertex_around_target_circulator vvc(m.halfedge(v),m), vvce(vvc); 
  BOOST_CHECK_EQUAL(CGAL::circulator_distance(vvc, vvce), 4);

  vvc =  vvce = Sm::Vertex_around_target_circulator(m.halfedge(u),m);
  BOOST_CHECK_EQUAL(CGAL::circulator_distance(vvc, vvce), 2);

  vvc = vvce = Sm::Vertex_around_target_circulator(m.halfedge(w),m);
  BOOST_CHECK_EQUAL(CGAL::circulator_distance(vvc, vvce), 3);

  vvc = vvce = Sm::Vertex_around_target_circulator(m.halfedge(x),m);
  BOOST_CHECK_EQUAL(CGAL::circulator_distance(vvc, vvce), 3);

  Sm::Halfedge_around_target_circulator hvc(m.halfedge(v),m), hvce(hvc); 
  BOOST_CHECK_EQUAL(CGAL::circulator_distance(hvc, hvce), 4);
  
  Sm::Face_around_target_circulator fvc(m.halfedge(v),m), fvce(fvc);
  BOOST_CHECK_EQUAL(CGAL::circulator_distance(fvc, fvce), 4);

  Sm::Vertex_around_face_circulator vfc(m.halfedge(f1),m), vfce(vfc);
  BOOST_CHECK_EQUAL(CGAL::circulator_distance(vfc, vfce), 3);

  Sm::Halfedge_around_face_circulator hfc(m.halfedge(f1),m), hfce(hfc); 
  BOOST_CHECK_EQUAL(CGAL::circulator_distance(hfc, hfce), 3);
}

BOOST_FIXTURE_TEST_CASE(emptiness, Surface_fixture)
{
  // add an isolated vertex
  Sm::Vertex_descriptor iv = m.add_vertex(Point_3(2,2,0));
  BOOST_CHECK(m.is_isolated(iv));

  Sm::Vertex_around_target_range vr = m.vertices_around_target(m.halfedge(iv));
  BOOST_CHECK(is_empty_range(boost::begin(vr), boost::end(vr)));

  Sm::Face_around_target_range fr = m.faces_around_target(m.halfedge(iv));
  BOOST_CHECK(is_empty_range(boost::begin(fr), boost::end(fr)));

  Sm::Halfedge_around_target_range hr = m.halfedges_around_target(m.halfedge(iv));
  BOOST_CHECK(is_empty_range(boost::begin(hr), boost::end(hr)));
  // not true for everything else
  m.remove_vertex(iv);
  BOOST_CHECK(m.is_removed(iv));
  Sm::Vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = m.vertices(); vb != ve; ++vb) {
    Sm::Vertex_around_target_range vr = m.vertices_around_target(m.halfedge(*vb));
    BOOST_CHECK(!is_empty_range(boost::begin(vr), boost::end(vr)));
  }
}

// trick cgal_create_CMakeLists
// int main()
// {
//   return 0;
// }
