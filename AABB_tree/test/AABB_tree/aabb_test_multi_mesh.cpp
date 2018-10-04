#include <fstream>
#include <algorithm>
#include <iterator>

#include <boost/functional/value_factory.hpp>
#include <boost/array.hpp>

#include <CGAL/assertions.h>
#include <CGAL/algorithm.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
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
CGAL::Tag_false> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Primitive_id Primitive_id;
typedef CGAL::Timer Timer;

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
  Tree tree(faces(m1).first, faces(m1).second, m1);
  tree.insert(faces(m2).first, faces(m2).second, m2);
  //tree.insert(faces(m1).first, faces(m1).second, m1);
  tree.build();
  Tree::Bounding_box bbox = tree.bbox();
  Point bbox_center((bbox.xmin() + bbox.xmax()) / 2, 
                     (bbox.ymin() + bbox.ymax()) / 2, 
                     (bbox.zmin() + bbox.zmax()) / 2);
  std::vector< Primitive_id > intersections;
  Ray ray(bbox_center+Vector(3,-0.25,0),bbox_center+Vector(-3,+0.25,0));
  tree.all_intersected_primitives(ray,
                                  std::back_inserter(intersections));
  
  return 0;
}
