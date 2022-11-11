
#include <iostream>
#include <fstream>
#include <CGAL/Timer.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Epick K;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh,
    CGAL::Default, CGAL::Tag_true, CGAL::Tag_true> Primitive2;
typedef CGAL::AABB_traits<K, Primitive2> Traits2;
typedef CGAL::AABB_tree<Traits2> Tree2;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh,
    CGAL::Default, CGAL::Tag_false, CGAL::Tag_true> Primitive3;
typedef CGAL::AABB_traits<K, Primitive3> Traits3;
typedef CGAL::AABB_tree<Traits3> Tree3;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh,
    CGAL::Default, CGAL::Tag_false, CGAL::Tag_false> Primitive4;
typedef CGAL::AABB_traits<K, Primitive4> Traits4;
typedef CGAL::AABB_tree<Traits4> Tree4;

int main(void)
{
  Mesh m;
  std::ifstream in("data/cube.off");
  if(in)
    in >> m;
  else{
    std::cout << "error reading bunny" << std::endl;
    return 1;
  }

  Tree t1(faces(m).begin(), faces(m).end(), m);
  Tree2 t2(faces(m).begin(), faces(m).end(), m);
  Tree3 t3(faces(m).begin(), faces(m).end(), m);
  Tree4 t4(faces(m).begin(), faces(m).end(), m);

  t1.build();
  t2.build();
  t3.build();
  t4.build();

  Primitive p1(faces(m).begin(), m);
  Primitive2 p2(faces(m).begin(), m);
  Primitive3 p3(faces(m).begin(), m);
  Primitive4 p4(faces(m).begin(), m);
  Triangle tr1 = t1.datum(p1);
  Triangle tr2 = t2.datum(p2);
  Triangle tr3 = t3.datum(p3);
  Triangle tr4 = t4.datum(p4);
  if(tr1 != tr2
     || tr1 != tr3
     || tr1 != tr4
     || tr2 != tr3
     || tr2 != tr4
     || tr3 != tr4)
    return 1;
  return 0;
}

