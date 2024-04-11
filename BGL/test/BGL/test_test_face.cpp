#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/generators.h>

#include <vector>
#include <iostream>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> SM;

typedef boost::graph_traits<SM>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<SM>::vertex_iterator vertex_iterator;
typedef std::vector<vertex_descriptor> V;
int main()
{
  {
    SM sm;

    vertex_descriptor vp = CGAL::add_vertex(sm);
    vertex_descriptor vq = CGAL::add_vertex(sm);
    vertex_descriptor vr = CGAL::add_vertex(sm);
    vertex_descriptor vs = CGAL::add_vertex(sm);
    std::array<vertex_descriptor,0> face0;
    assert( ! CGAL::Euler::can_add_face(face0,sm) );
    std::array<vertex_descriptor,1> face1 = { vp };
    assert( ! CGAL::Euler::can_add_face(face1,sm) );
    std::array<vertex_descriptor,2> face2 = { vp, vq };
    assert( ! CGAL::Euler::can_add_face(face2,sm) );

    std::array<vertex_descriptor,3> face = { vp, vq, vr };
    CGAL::Euler::add_face(face, sm);

    assert( ! CGAL::Euler::can_add_face(face,sm) );
    std::swap(face[0],face[1]);
    assert( CGAL::Euler::can_add_face(face,sm) );

    face[2] = vs;
    assert( CGAL::Euler::can_add_face(face,sm) );
    std::swap(face[0],face[1]);
    assert( ! CGAL::Euler::can_add_face(face,sm) );
  }

  {
    SM sm;
    Point_3 p(0,0,0), q(1,0,0), r(0,1,0), s(0,0,1);
    CGAL::make_tetrahedron(p, q, r, s, sm);
    std::array<vertex_descriptor,3> face;
    vertex_iterator it = vertices(sm).first;
    face[0] = *it;
    ++it;
    face[1] = *it;
    face[2] = CGAL::add_vertex(sm);
    assert( ! CGAL::Euler::can_add_face(face,sm) );
    std::swap(face[0],face[1]);
    assert( ! CGAL::Euler::can_add_face(face,sm) );
  }

  std::cout << "Done" << std::endl;

  return 0;
}
