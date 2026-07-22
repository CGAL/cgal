#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Ray_3 Ray;

typedef CGAL::Surface_mesh<Point> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits_3<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef std::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;

struct Skip
{
  face_descriptor fd;

  Skip(const face_descriptor fd)
    : fd(fd)
  {}

  bool operator()(const face_descriptor& t) const
  { if(t == fd){
      std::cerr << "ignore " << t  <<std::endl;
    };
    return(t == fd);
  }

};

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/tetrahedron.off");

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  Tree tree(faces(mesh).first, faces(mesh).second, mesh);

  double d = CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)?-1:1;

  for(face_descriptor fd : faces(mesh))
  {
    halfedge_descriptor hd = halfedge(fd,mesh);
    Point p = CGAL::centroid(mesh.point(source(hd,mesh)),
                             mesh.point(target(hd,mesh)),
                             mesh.point(target(next(hd,mesh),mesh)));
    Vector v = CGAL::Polygon_mesh_processing::compute_face_normal(fd,mesh);

    Ray ray(p,d * v);
    Skip skip(fd);
    Ray_intersection intersection = tree.first_intersection(ray, skip);
    if(intersection)
    {
      if(std::get_if<Point>(&(intersection->first))){
        const Point* p =  std::get_if<Point>(&(intersection->first) );
        std::cout <<  *p << std::endl;
      }
    }
  }

  std::cerr << "done" << std::endl;

  return 0;
}
