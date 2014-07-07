#define BOOST_TEST_MODULE graph_geometry test
#include <boost/test/unit_test.hpp>

#include "test_Prefix.h"

#include <CGAL/boost/graph/Graph_geometry.h>

// BOOST_AUTO_TEST_CASE( compute_face_normals_test )
// {
//   Surface_fixture f;
//   Cube_fixture cf;
//   f.m.add_property<SM::Face, K::Vector_3>("f:normal");
//   CGAL::Property<SM::Face, K::Vector_3> 
//     facenormals = f.m.get_property<SM::Face, K::Vector_3>("f:normal");

//   CGAL::calculate_face_normals(f.m, boost::get(CGAL::vertex_point, f.m), facenormals);
// }

// BOOST_AUTO_TEST_CASE( compute_vertex_normals_test )
// {
//   Surface_fixture f;
//   f.m.add_property<SM::Vertex, K::Vector_3>("v:normal");

//   CGAL::Property<SM::Vertex, K::Vector_3> 
//     normals = f.m.get_property<SM::Vertex, K::Vector_3>("v:normal");

//   CGAL::calculate_vertex_normals(f.m, boost::get(CGAL::vertex_point, f.m), 
//                                  normals, boost::get(CGAL::vertex_is_border, f.m));
// }

BOOST_AUTO_TEST_CASE( triangle_test )
{
  Polyhedron p;
  {
    std::ifstream in("data/7_faces_triangle.off");
    BOOST_CHECK(in >> p);
  }
  typedef boost::graph_traits<Polyhedron> Traits;
  typedef Traits::face_iterator face_iterator;
  typedef Traits::face_iterator edge_descriptor;
  typedef Traits::vertex_descriptor vertex_descriptor;

  face_iterator fb, fe;
  for(boost::tie(fb, fe) = faces(p); fb != fe; ++fb) {
    // call the overloaded function
    Kernel::Triangle_3 tri = CGAL::triangle(p, halfedge(*fb, p));
      // call the actual impl
    Kernel::Triangle_3 tri2 = CGAL::triangle(p, halfedge(*fb, p), get(CGAL::vertex_point, p));
    // fiddle up a triangle ourselves
    std::vector<vertex_descriptor> vs;
    vs.push_back(target(halfedge(*fb, p),p));
    vs.push_back(target(next(halfedge(*fb, p),p),p));
    vs.push_back(target(next(next(halfedge(*fb, p),p),p),p));

    BOOST_CHECK_EQUAL(vs.size(), 3);
    Kernel::Triangle_3 tricomp = Kernel::Triangle_3(get(CGAL::vertex_point, p, vs[0]), 
                                                    get(CGAL::vertex_point, p, vs[1]),
                                                    get(CGAL::vertex_point, p, vs[2]));
    BOOST_CHECK_EQUAL(tri, tricomp);
    BOOST_CHECK_EQUAL(tri2, tricomp);
    BOOST_CHECK_EQUAL(tri, tri2);
  }
}

BOOST_AUTO_TEST_CASE( is_pure_triangle_test )
{
  Polyhedron pure, impure;
  {
    std::ifstream in("data/7_faces_triangle.off");
    BOOST_CHECK(in >> pure);
  }
  {
    // we have no impure data, this is a major snafu
  }
  BOOST_CHECK(CGAL::is_pure_triangle(pure));
}

// int main()
// {
//   return 0;
// }
