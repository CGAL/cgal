
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
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh,
    CGAL::Default, CGAL::Tag_true, CGAL::Tag_true> Primitive_cached;

typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_traits<K, Primitive_cached> Traits_cached;
typedef CGAL::AABB_tree<Traits> Tree;
typedef CGAL::AABB_tree<Traits_cached> Tree_cached;
typedef Tree::Primitive_id Primitive_id;
typedef Tree_cached::Primitive_id Primitive_id_cached;

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
  Tree_cached t2(faces(m).begin(), faces(m).end(), m);
  
  t1.build();
  t2.build();
  Primitive p1(faces(m).begin(), m);
  Primitive_cached p2(faces(m).begin(), m);
  Triangle tr1 = t1.datum(p1);
  Triangle tr2 = t2.datum(p2);
  if(tr1 != tr2)
    return 1;
  return 0;
}

