#include <fstream>
#include <algorithm>
#include <iterator>

#include <boost/functional/value_factory.hpp>

#include <CGAL/algorithm.h>
#include <CGAL/point_generators_3.h>
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
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh,
CGAL::Default,
CGAL::Tag_false> T_Primitive;
typedef CGAL::AABB_traits<K, T_Primitive> T_Traits;
typedef CGAL::AABB_tree<T_Traits> T_Tree;
typedef T_Tree::Primitive_id T_Primitive_id;

typedef CGAL::AABB_halfedge_graph_segment_primitive<Mesh,
CGAL::Default,
CGAL::Tag_false> E_Primitive;
typedef CGAL::AABB_traits<K, E_Primitive> E_Traits;
typedef CGAL::AABB_tree<E_Traits> E_Tree;
typedef E_Tree::Primitive_id E_Primitive_id;

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
  in.open(CGAL::data_file_path("meshes/tetrahedron.off"));
  if(in)
    in >> m2;
  else{
    std::cout << "error reading tetrahedron" << std::endl;
    return 1;
  }
  in.close();
  T_Tree tree(faces(m1).first, faces(m1).second, m1);
  tree.insert(faces(m2).first, faces(m2).second, m2);
  tree.build();
  T_Tree::Bounding_box bbox = tree.bbox();
  Point bbox_center((bbox.xmin() + bbox.xmax()) / 2,
                     (bbox.ymin() + bbox.ymax()) / 2,
                     (bbox.zmin() + bbox.zmax()) / 2);
  std::vector< T_Primitive_id > intersections;
  Ray ray(bbox_center+Vector(3,-0.25,0),bbox_center+Vector(-3,+0.25,0));
  tree.all_intersected_primitives(ray,
                                  std::back_inserter(intersections));
  E_Tree e_tree(edges(m1).first, edges(m1).second, m1);
  e_tree.insert(edges(m2).first, edges(m2).second, m2);
  e_tree.build();
  std::vector< E_Primitive_id > e_intersections;
  Ray e_ray(Point(0,0,0),Point(0,1,1));
  e_tree.all_intersected_primitives(e_ray,
                                  std::back_inserter(e_intersections));



  return 0;
}
