#define BOOST_TEST_MODULE surface_mesh test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include "SM_common.h"

#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>

#include <boost/bind.hpp>
#include <boost/range/algorithm.hpp>

BOOST_AUTO_TEST_CASE( constructors_test )
{
  Sm mesh;
  Sm mesh2(mesh);
  mesh = mesh2;
  // Surface_mesh::assign
}

template<typename Iterator>
void test_iterator(Iterator begin, Iterator end, 
                   typename std::iterator_traits<Iterator>::difference_type distance) 
{
  BOOST_CHECK(begin != end);
  BOOST_CHECK(begin == begin);
  BOOST_CHECK(std::distance(begin, end) == distance);
}

BOOST_AUTO_TEST_CASE( standard_iterators ) 
{
  Surface_fixture f;

  Sm::Vertex_iterator vb, ve;
  boost::tie(vb, ve) = f.m.vertices();
  test_iterator(vb, ve, 5);

  Sm::Halfedge_iterator hb, he;
  boost::tie(hb, he) = f.m.halfedges();
  test_iterator(hb, he, 14);

  Sm::Edge_iterator eb, ee;
  boost::tie(eb, ee) = f.m.edges();
  test_iterator(eb, ee, 7);

  Sm::Face_iterator fb, fe;
  boost::tie(fb, fe) = f.m.faces();
  test_iterator(fb, fe, 3);
}


BOOST_AUTO_TEST_CASE( const_iterators ) {
  Surface_fixture f;
}

BOOST_AUTO_TEST_CASE( iter_inteop ) { 
  /* unimplemented */
}

BOOST_AUTO_TEST_CASE( test_descriptors ) {
  Sm::Vertex_index v;
  Sm::Halfedge_index h;
  v < v;
  h == h;
  // does not compile and so it should
  // v < h;
}

BOOST_AUTO_TEST_CASE( test_remove_edge ) 
{
  Surface_fixture f;
  Sm::Halfedge_index wv, uw, wx, xv, vu;
  wv = f.m.halfedge(f.w, f.v);
  BOOST_CHECK_EQUAL(f.m.target(wv), f.v);
  BOOST_CHECK(wv.is_valid());
  uw = f.m.halfedge(f.u, f.w);
  BOOST_CHECK(uw.is_valid());
  wx = f.m.halfedge(f.w, f.x);
  BOOST_CHECK(wx.is_valid());
  xv = f.m.halfedge(f.x, f.v);
  BOOST_CHECK(xv.is_valid());
  vu = f.m.halfedge(f.v, f.u);
  BOOST_CHECK(vu.is_valid());
  BOOST_CHECK_EQUAL(f.m.degree(f.v), 4);
  BOOST_CHECK_EQUAL(f.m.degree(f.w), 3);

  f.m.set_next(uw, wx);
  f.m.set_next(xv, vu);
  f.m.remove_edge(Sm::Edge_index(wv));
  Sm::Halfedge_around_target_circulator a(f.m.halfedge(f.w),f.m), b(a);
  BOOST_CHECK_EQUAL(CGAL::circulator_distance(a, b), 2);
  BOOST_CHECK_EQUAL(f.m.degree(f.w), 2);
  a = b = Sm::Halfedge_around_target_circulator(f.m.halfedge(f.v),f.m);
  BOOST_CHECK_EQUAL(CGAL::circulator_distance(a, b), 3);
  BOOST_CHECK_EQUAL(f.m.degree(f.v), 3);

  // now remove a border edge to check if this works
  // this should not lower the number of faces
  Sm::size_type old_removed_faces_size = f.m.number_of_removed_faces();
  f.m.remove_edge(Sm::Edge_index(wx));
  BOOST_CHECK_EQUAL(f.m.number_of_removed_faces(), old_removed_faces_size);
}


BOOST_AUTO_TEST_CASE( memory_reuse_test )
{
  Cube_fixture f;
  // buffer all faces
  typedef std::vector<Sm::Vertex_index> VecFace;
  typedef std::vector<VecFace> Faces;

  Faces faces;
  Sm::Face_iterator fb, fe;
  for(boost::tie(fb, fe) = f.m.faces(); fb != fe; ++fb) {
    faces.push_back(VecFace());
    Sm::Vertex_around_face_circulator vafb(f.m.halfedge(*fb), f.m), vafe(vafb);
    if(vafb)
      do {
        faces.back().push_back(*vafb);
      } while(++vafb != vafe);
    BOOST_CHECK_EQUAL(f.m.degree(*fb), faces.back().size());
  }

  Sm::Vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = f.m.vertices(); vb != ve; ++vb) {
    f.m.set_halfedge(*vb, Sm::Halfedge_index());
  }
  
  // remove all faces
  std::size_t old_face_size = f.m.number_of_faces(); 
  std::size_t old_removed_face_size = f.m.number_of_removed_faces();
  boost::range::for_each(f.m.faces(), boost::bind(&Sm::remove_face, boost::ref(f.m), _1));
  BOOST_CHECK_EQUAL(f.m.number_of_faces(), 0);
  BOOST_CHECK_EQUAL(f.m.number_of_removed_faces(), old_face_size + old_removed_face_size);
  // remove all edges
  std::size_t old_edge_size = f.m.number_of_edges();
  std::size_t old_removed_edge_size = f.m.number_of_removed_edges();
  boost::range::for_each(f.m.edges(), 
                         boost::bind(static_cast<void (Sm::*)(Sm::Edge_index)>(&Sm::remove_edge), 
                                     boost::ref(f.m), _1));
  BOOST_CHECK_EQUAL(f.m.number_of_faces() , 0);
  BOOST_CHECK_EQUAL(f.m.number_of_removed_edges(), old_edge_size + old_removed_edge_size);

  int fc = 0;
  // add all again
  for(Faces::iterator it = faces.begin(); it != faces.end(); ++it) {
    Sm::Face_index fd = f.m.add_face(*it);
    BOOST_CHECK(fd.is_valid());
    f.m.set_vertex_halfedge_to_border_halfedge(f.m.halfedge(fd));
    for(VecFace::iterator it2 = it->begin(); it2 != it->end(); ++it2) { 

      Sm::Halfedge_index h = f.m.halfedge(*it2);
      Sm::Face_index fa = f.m.face(h);

    }
  }
  
  BOOST_CHECK_EQUAL(f.m.number_of_edges() , old_edge_size);
  BOOST_CHECK_EQUAL(f.m.number_of_faces() , old_face_size);


  // remove all vertices
  std::size_t old_size = f.m.number_of_vertices();
  std::size_t old_removed_size = f.m.number_of_removed_vertices();

  boost::range::for_each(f.m.vertices(), boost::bind(&Sm::remove_vertex, boost::ref(f.m), _1));
  BOOST_CHECK_EQUAL(f.m.number_of_vertices() , 0);
  BOOST_CHECK_EQUAL(f.m.number_of_removed_vertices(), old_size + old_removed_size);

  for(unsigned int i = 0; i < old_size; ++i)
  {
    f.m.add_vertex(K::Point_3());
  }

  // the size must remain the same
  BOOST_CHECK_EQUAL(f.m.number_of_vertices(), old_size);
}

BOOST_AUTO_TEST_CASE( test_validate )
{
  Cube_fixture cf;
  Surface_fixture f1;
  Surface_fixture_2 f2;
  Surface_fixture_3 f3;
  BOOST_CHECK(cf.m.is_valid());
  BOOST_CHECK(f1.m.is_valid());
  BOOST_CHECK(f2.m.is_valid());  
  BOOST_CHECK(f3.m.is_valid());
}

BOOST_AUTO_TEST_CASE(isolated_vertex_check)
{
  Surface_fixture f;
  Sm::Vertex_index isolated = f.m.add_vertex(Point_3(10, 10, 10));
  BOOST_CHECK(f.m.is_isolated(isolated));
  BOOST_CHECK(!f.m.halfedge(isolated).is_valid());
  BOOST_CHECK(! f.m.is_border(isolated));
  BOOST_CHECK(f.m.degree(isolated) == 0);
}

BOOST_AUTO_TEST_CASE(embedded_vertex_check)
{
  Surface_fixture_2 f;
  BOOST_CHECK(!f.m.is_isolated(f.y));
  BOOST_CHECK(!f.m.is_border(f.y));
  BOOST_CHECK(f.m.halfedge(f.y).is_valid());
  BOOST_CHECK(f.m.degree(f.y) == 4);
}

BOOST_AUTO_TEST_CASE(border_vertex_check)
{
  Surface_fixture f;
  BOOST_CHECK(!f.m.is_isolated(f.y));
  BOOST_CHECK(f.m.is_border(f.y));
  BOOST_CHECK(f.m.degree(f.y) == 2);
  BOOST_CHECK(f.m.halfedge(f.y).is_valid());
}


BOOST_AUTO_TEST_CASE( point_position_accessor )
{
  Surface_fixture f;
  // by property
  f.m.points()[f.x];
  // non-const
  f.m.point(f.x);
  // const as an lvalue
  f.m.point(f.x) = CGAL::ORIGIN;
  BOOST_CHECK(f.m.point(f.x) == CGAL::ORIGIN);
}

BOOST_AUTO_TEST_CASE( properties ) {
  Surface_fixture f;
  f.m.add_property_map<Sm::Vertex_index, int>("illuminatiproperty", 23);
  // CGAL::Property_map<Sm::Vertex_index, int> prop = f.m.get_property_map<Sm::Vertex_index, int>("illuminatiproperty");
}


// trick cgal_create_CMakeLists
// int main()
// {
//   return 0;
// }
