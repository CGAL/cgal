#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point_3;
typedef CGAL::Bbox_3 Bbox_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;
typedef CGAL::Timer Timer;


template <typename TriangleMesh>
void triangle_mesh(const char* fname)
{
  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh> Primitive;
  typedef CGAL::AABB_traits<K, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;

  TriangleMesh tmesh;
  std::ifstream in(fname);
  in >> tmesh;
  Timer t;
  t.start();
  Tree tree(faces(tmesh).first, faces(tmesh).second, tmesh);
  tree.build();
  std::cout << t.time() << " sec." << std::endl;
  std::cout << "Closest point to ORIGIN:" << tree.closest_point(CGAL::ORIGIN) << std::endl;
}


Bbox_3 bbox(boost::graph_traits<Surface_mesh>::face_descriptor fd,
            const Surface_mesh& p)
{
  boost::graph_traits<Surface_mesh>::halfedge_descriptor hd = halfedge(fd,p);
  Bbox_3 res = p.point(source(hd,p)).bbox();
  res += p.point(target(hd,p)).bbox();
  res += p.point(target(next(hd,p),p)).bbox();
  return res;
}


void surface_mesh_cache_bbox(const char* fname)
{
  typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;
  typedef Surface_mesh::Property_map<face_descriptor,Bbox_3> Bbox_pmap;
  typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> Primitive;
  typedef CGAL::AABB_traits<K, Primitive,Bbox_pmap> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;

  Surface_mesh tmesh;
  std::ifstream in(fname);
  in >> tmesh;

  Timer t;
  t.start();
  Bbox_pmap bb = tmesh.add_property_map<face_descriptor,Bbox_3>("f:bbox",Bbox_3()).first;

  for(face_descriptor fd : faces(tmesh)){
    put(bb, fd, bbox(fd,tmesh));
  }
  Traits traits(bb);
  Tree tree(traits);
  tree.insert(faces(tmesh).first, faces(tmesh).second, tmesh);
  tree.build();
  tmesh.remove_property_map(bb);

  std::cout << t.time() << " sec."<< std::endl;
   std::cout << "Closest point to ORIGIN:" << tree.closest_point(CGAL::ORIGIN) << std::endl;
}


int main(int argc, char* argv[])
{
  std::cout << "Polyhedron_3" << std::endl;
  triangle_mesh<Polyhedron_3>((argc>1)?argv[1]:"data/tetrahedron.off");

  std::cout << "Surface_mesh" << std::endl;
  triangle_mesh<Surface_mesh>((argc>1)?argv[1]:"data/tetrahedron.off");

  std::cout << "Surface_mesh with cached Bbox_3" << std::endl;
  surface_mesh_cache_bbox((argc>1)?argv[1]:"data/tetrahedron.off");

  return EXIT_SUCCESS;
}
