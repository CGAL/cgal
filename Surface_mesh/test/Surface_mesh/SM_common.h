#ifndef CGAL_SM_COMMON_H
#define CGAL_SM_COMMON_H

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO.h>

#include <CGAL/Simple_cartesian.h>

#include <boost/assign.hpp>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Sm;
typedef K::Point_3 Point_3;

/* Mesh used in this fixture:
//
// yo
//  |\
//  | \
//  |f3\
// vo---ox
//  |\f2|
//  | \ |
//  |f1\|
// uo---ow
*/
struct Surface_fixture {
  Surface_fixture() {
    u = m.add_vertex(Point_3(0,0,0));
    v = m.add_vertex(Point_3(1,0,0));
    w = m.add_vertex(Point_3(0,1,0));
    x = m.add_vertex(Point_3(1,1,0));
    y = m.add_vertex(Point_3(2,0,0));
    f1 = m.add_face(u,w,v);
    f2 = m.add_face(v,w,x);
    f3 = m.add_face(v,x,y);
  }

  Sm m;
  Sm::Vertex_index u, v, w, x, y;
  Sm::Face_index f1, f2, f3;

  ~Surface_fixture() {}
};


struct Surface_fixture_2 {

  // u        v
  // +--------+
  // |\      /|
  // | \ f2 / |
  // |  \y /  |
  // | f3\/ f1|
  // |   /\   |
  // |  /  \  |
  // | / f4 \ |
  // |/      \|
  // +--------+
  // w        x
  Surface_fixture_2() {
    u = m.add_vertex(Point_3(0,2,0));
    v = m.add_vertex(Point_3(2,2,0));
    w = m.add_vertex(Point_3(0,0,0));
    x = m.add_vertex(Point_3(2,0,0));
    y = m.add_vertex(Point_3(1,1,0));
    f1 = m.add_face(x, v, y);
    f2 = m.add_face(u, y, v);
    f3 = m.add_face(u, w, y);
    f4 = m.add_face(w, x, y);
  }

  Sm m;
  Sm::Vertex_index u, v, w, x, y;
  Sm::Face_index f1, f2, f3, f4;

  ~Surface_fixture_2() {}
};

struct Surface_fixture_3 {

   // u            x            z
   // +------------+------------+
   // |            |            |
   // |            |            |
   // |      f1    |    f2      |
   // |            |            |
   // |            |            |
   // +------------+------------+
   // v            w            y
  Surface_fixture_3() {
    u = m.add_vertex(Point_3(0,1,0));
    v = m.add_vertex(Point_3(0,0,0));
    w = m.add_vertex(Point_3(1,0,0));
    x = m.add_vertex(Point_3(1,1,0));
    y = m.add_vertex(Point_3(2,0,0));
    z = m.add_vertex(Point_3(2,1,0));

    std::vector<Sm::Vertex_index> vec;
    using namespace boost::assign;
    vec += u, v, w, x;
    f1 = m.add_face(vec);
    vec.clear();
    vec += x, w, y, z;
    f2 = m.add_face(vec);
  }

  Sm m;
  Sm::Vertex_index u, v, w, x, y, z;
  Sm::Face_index f1, f2;

  ~Surface_fixture_3() {}
};


struct Cube_fixture {
  Cube_fixture() {
    CGAL::read_mesh(m, "cube.off");
  }

  Sm m;

  ~Cube_fixture() {}
};


#endif /* CGAL_SM_COMMON_H */
