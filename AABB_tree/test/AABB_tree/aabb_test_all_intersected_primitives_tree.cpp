#include <fstream>
#include <iterator>

#include <CGAL/assertions.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/Timer.h>

typedef CGAL::Epick K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Segment_3 Segment;
typedef K::Ray_3 Ray;
typedef CGAL::Surface_mesh<CGAL::Point_3<CGAL::Epick> > Mesh;
typedef CGAL::AABB_halfedge_graph_segment_primitive<Mesh,
CGAL::Default,
CGAL::Tag_false> S_Primitive;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh,
CGAL::Default,
CGAL::Tag_false> T_Primitive;
typedef CGAL::AABB_traits<K, T_Primitive> T_Traits;
typedef CGAL::AABB_traits<K, S_Primitive> S_Traits;
typedef CGAL::AABB_tree<T_Traits> T_Tree;
typedef CGAL::AABB_tree<S_Traits> S_Tree;
typedef T_Tree::Primitive_id T_Primitive_id;
typedef S_Tree::Primitive_id S_Primitive_id;

int main()
{
  CGAL::Surface_mesh<CGAL::Point_3<CGAL::Epick> > m1, m2;
  std::ifstream in("data/cube.off");
  if(in)
    in >> m1;
  else{
    std::cout << "error reading cube" << std::endl;
    return 1;
  }
  in.close();
  in.open("data/tetrahedron.off");
  if(in)
    in >> m2;
  else{
    std::cout << "error reading tetrahedron" << std::endl;
    return 1;
  }
  in.close();
  T_Tree cube_tree(faces(m1).first, faces(m1).second, m1);
  S_Tree tet_tree(edges(m2).first, edges(m2).second, m2);
  cube_tree.build();
  tet_tree.build();

  std::list<T_Tree::Primitive::Id> t_primitives;
  std::list<S_Tree::Primitive::Id> s_primitives;
  cube_tree.all_intersected_primitives(tet_tree,std::back_inserter(t_primitives));
  CGAL_assertion(t_primitives.size() == 6);
  tet_tree.all_intersected_primitives(cube_tree,std::back_inserter(s_primitives));
  CGAL_assertion(s_primitives.size() == 6);
  CGAL_assertion(tet_tree.do_intersect(cube_tree));
  CGAL_assertion(cube_tree.do_intersect(tet_tree));

  std::vector<T_Tree::Primitive::Id> all_primitives;
  cube_tree.all_intersected_primitives(tet_tree, std::back_inserter(all_primitives));
  bool found_f5 = false;
  for(auto prim : all_primitives)
  {
    if((int)prim.first == 5)
      found_f5 = true;
  }
  CGAL_assertion(found_f5);
  CGAL_USE(found_f5);
  return 0;
}
