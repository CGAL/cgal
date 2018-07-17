#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_transformed_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/centroid.h>


typedef CGAL::Epick K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Segment_3 Segment;
typedef K::Ray_3 Ray;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Base_traits;
typedef CGAL::AABB_transformed_traits<Base_traits, K> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Primitive_id Primitive_id;

int main()
{
  typedef std::vector< Tree::Intersection_and_primitive_id<Ray>::Type > IntersectionVector;
  Polyhedron polyhedron;
  std::ifstream in("data/bunny00.off");
  if(in)
    in >> polyhedron;
  else{
    std::cout << "error reading bunny" << std::endl;
    return 1;
  }
  CGAL::Aff_transformation_3<K> transfo(CGAL::TRANSLATION, Vector(250,300,1));
  IntersectionVector all_intersections;
  Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
  tree.traits().set_transformation(transfo);
  tree.build();
  tree.traits().set_transformation(CGAL::Aff_transformation_3<K>(CGAL::IDENTITY));
  Point bunny_center =
      CGAL::centroid(polyhedron.points_begin(), polyhedron.points_end());
      
  Ray ray1(Point(0,0,-10),bunny_center);
  tree.all_intersections(ray1, std::back_inserter(all_intersections));
  std::size_t nb_inter = all_intersections.size();
  std::cout<<nb_inter<<" primitives intersected before translation."<<std::endl;
  all_intersections.clear();
  tree.traits().set_transformation(transfo);
  CGAL::Polygon_mesh_processing::transform(transfo,polyhedron);
  
  Ray ray2(transfo.transform(Point(0,0,-10)),transfo.transform(bunny_center));
  tree.all_intersections(ray2, std::back_inserter(all_intersections));
  std::cout<<all_intersections.size()<<" primitives intersected after translation."<<std::endl;
  if(nb_inter == all_intersections.size())
    return 0;
  
  return 1;
}
