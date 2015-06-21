
#include "test_Prefix.h"

#include <CGAL/boost/graph/Graph_geometry.h>


void triangle_test()
{
  Polyhedron p;
  {
    std::ifstream in("data/7_faces_triangle.off");
    assert(in >> p);
  }
  typedef boost::graph_traits<Polyhedron> Traits;
  typedef Traits::face_iterator face_iterator;
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

    assert(vs.size() == 3);
    Kernel::Triangle_3 tricomp = Kernel::Triangle_3(get(CGAL::vertex_point, p, vs[0]), 
                                                    get(CGAL::vertex_point, p, vs[1]),
                                                    get(CGAL::vertex_point, p, vs[2]));
    assert(tri == tricomp);
    assert(tri2 == tricomp);
    assert(tri == tri2);
  }
}

void
is_pure_triangle_test()
{
  Polyhedron pure, impure;
  {
    std::ifstream in("data/7_faces_triangle.off");
    assert(in >> pure);
  }
  {
    // we have no impure data, this is a major snafu
  }
  assert(CGAL::is_triangle_mesh(pure));
}

int main()
{
  triangle_test();
  is_pure_triangle_test();
  
  std::cerr << "done\n";
  return 0;
}
